# Rui Shao, 2022 Sep
# Baum-Welch algorithm for discrete Hidden Markov models

# -----------------------------------------------------------------------------

# genomic region segmentation

# -----------------------------------------------------------------------------

#' Baum-welch HMM with log-normal probility estimation
#' @param observation A matrix of features in normal distribution
#' @param num_states  Number of hidden states
#' @param init_method Initiation HMM states with kmeans, mclust, or quantile methods.
#' 
run_hmm <- function(observation, num_states = 10, init_method = "kmeans") {
  
  # -------------------------------- Steps ------------------------------------
  #
  # 0. initializer()        # random log-normal centroids
  #
  # 1. baum_welch() 
  #   |
  #   |---- forward() and backward() to get alpha and beta
  #   |
  #   `---- gamma and xi to EM update transition_prob and emission_prob
  #
  # 2. viterbi() to get the optimal state path
  #
  # ---------------------------------------------------------------------------

  get_emission_probs <- function(observation, emission_paras) {
    
    # output: (state, position) matrix in log
    
    num_states = nrow(emission_paras[["mu"]])
    emission_probs <- matrix(1, nrow = num_states, ncol = nrow(observation))
    
    for (i in seq_len(num_states)) { # joint probability
      emission_probs[i, ] = 
        colSums(log(dnorm(t(observation), 
                          mean = emission_paras[["mu"]][i, ], 
                          sd = emission_paras[["sd"]][i, ]) + 1e-10))
    }
    rownames(emission_probs) <- rownames(emission_paras[["mu"]])
    emission_probs
  }
  
  
  initializer <- function(num_states = 10, observation, 
                          init_method = c("kmeans", "mclust", "quantile")) {
    # Args:
    # num_states: number of hidden states
    # observation: input observation matrix, log-transformed
    
    num_features = ncol(observation)
    states = seq_len(num_states)
    
    # use clustering method to initialize state means
    non_zero_observation = observation[rowSums(observation) > 0, ]
    if (any(init_method == "kmeans")) {
      message("Initialize HMM states with kmeans.")
      k_res = suppressWarnings(kmeans(non_zero_observation, 
                                      centers = num_states - 1,
                                      iter.max = 100))
      res_mu = k_res$centers
    } else if (any(init_method == "mclust")) {
      require(mclust)
      message("Initialize HMM states with multi-normal distribution.")
      m_res = mclust::Mclust(non_zero_observation, G = num_states - 1)
      res_mu = t(m_res$parameters$mean)
    } else {
      res_mu = apply(observation, 2,
                               function(x)
                                 sample(quantile(x, seq(0.1, 1, length.out = num_states - 1)))
                     )
    }
    
    # set initialize state means
    initi_state_mean = matrix(0, nrow = num_states, ncol = num_features)
    initi_state_mean[-num_states, ] = res_mu
    
    rownames(initi_state_mean) = states
    if (!is.null(colnames(observation))) colnames(initi_state_mean) = colnames(observation)
    
    # set initial state sd as the total
    initi_state_sd = apply(observation, 2, 
                           function(x) 
                             rep(sd(x), num_states) 
    ) 
    rownames(initi_state_sd) = rownames(initi_state_mean)
    colnames(initi_state_sd) = colnames(initi_state_mean)
    
    # states X observs
    emission_paras = list("mu" = initi_state_mean,  # mean of each state
                          "sd" = initi_state_sd)  # standard deviation of each state
    log_emission_probs = get_emission_probs(observation, emission_paras) # log probs
    
    transition_probs = matrix(1, nrow = num_states, ncol = num_states) / num_states
    dimnames(transition_probs) = list(states, states)
    
    # create a parameter list 
    hmm_param = list()
    hmm_param$states = states
    hmm_param$num_features = num_features
    
    hmm_param$initstate_probs = rep(1, num_states) / num_states
    hmm_param$log_transition_probs = log(transition_probs)
    
    hmm_param$emission_paras = emission_paras
    hmm_param$log_emission_probs = log_emission_probs
    hmm_param
  }
  
  
  forward <- function(observation, hmm_param) {
    # Args:
    # observation: position X feature matrix
    # hmm_param: a parameter list
    
    states = hmm_param$states
    log_initstate_probs  = log(hmm_param$initstate_probs)
    log_transition_probs = hmm_param$log_transition_probs
    log_emission_probs = hmm_param$log_emission_probs
    
    .len = nrow(observation)
    
    alpha = matrix(0, ncol = length(states), nrow = .len)
    colnames(alpha) = states
    
    # output in log
    alpha[1, ] = log_initstate_probs + log_emission_probs[, 1]
    
    for (i in seq_len(.len)[-1]) {
      for (j in states) { # 5 times slower
        
        alpha[i, j] = matrixStats::logSumExp(alpha[i-1, ] + 
                                               log_transition_probs[, j] + 
                                               log_emission_probs[j, i])
      }
    }
    alpha
  }
  
  
  backward <- function(observation, hmm_param) {
    # Args:
    # observation: position X feature matrix
    # hmm_param: a parameter list
    
    states =  hmm_param$states
    log_emission_probs   = hmm_param$log_emission_probs
    log_transition_probs = hmm_param$log_transition_probs
    
    .len = nrow(observation)
    beta = matrix(0, ncol = length(states), nrow = .len)
    colnames(beta) = states
    
    # output in log
    for (i in seq(.len - 1, 1)) {
      for (j in states) { 
        beta[i, j] = matrixStats::logSumExp(beta[i + 1, ] +
                                              log_transition_probs[j, ] + 
                                              log_emission_probs[, i + 1])
      }
    }
    beta
  }
  
  
  viterbi <- function(observation, hmm_param) {
    # Args:
    # observation: position X feature matrix
    # hmm_param: a parameter list
    
    message("Find path with viterbi.")
    
    states =  hmm_param$states
    log_initstate_probs  = log(hmm_param$initstate_probs)
    log_emission_probs   = hmm_param$log_emission_probs
    log_transition_probs = hmm_param$log_transition_probs
    
    .len = nrow(observation)
    
    viterbi = `names<-`(numeric(length(states)), states)
    state_path = lapply(states, function(x) c(x, numeric(.len - 1)))
    names(state_path) = states
    
    viterbi = log_initstate_probs + log_emission_probs[, 1]
    
    for (t in seq_len(.len)[-1]) {
      
      tmp_path = state_path
      viterbi_tmp = `names<-`(numeric(length(states)), states)
      viterbi_prod = viterbi + log_transition_probs
      
      for (state_to in states) {
        intermediate_probs = viterbi_prod[, state_to]
        
        max_state = names(intermediate_probs)[which.max(intermediate_probs)]
        max_prob = max(intermediate_probs)
        
        prob = log_emission_probs[state_to, t] + max_prob   
        viterbi_tmp[state_to] = prob
        
        tmp_state_path = state_path[[max_state]]
        tmp_state_path[t] = state_to
        tmp_path[[state_to]] = tmp_state_path
      }
      
      viterbi = viterbi_tmp
      state_path = tmp_path
    }
    
    max_path = state_path[[which.max(viterbi)]]
    
    list(max(viterbi), max_path)
  }
  
  
  baum_welch <- function(observation, hmm_param) {
    # Args:
    # observation: position X feature matrix
    # hmm_param: a parameter list
    
    states =  hmm_param$states
    emission_paras = hmm_param$emission_paras
    log_emission_probs   = hmm_param$log_emission_probs
    log_transition_probs = hmm_param$log_transition_probs
    
    .len = nrow(observation)
    if (.len == 0) return(NULL)
    
    alpha = forward(observation, hmm_param)
    beta = backward(observation, hmm_param)  
    
    gamma = matrix(0, ncol = length(states), nrow = .len)
    xi = array(0, dim = c(.len - 1, length(states), length(states)), 
               dimnames = list(seq_len(.len - 1), states, states))
    colnames(gamma) = states
    
    # E step
    gamma = alpha + beta
    gamma = gamma - matrixStats::rowLogSumExps(gamma)
    
    for (t in seq_len(.len - 1)) {
      xi[t,,] = t(t(alpha[t, ] + 
                      log_transition_probs) + 
                    log_emission_probs[, t + 1] +
                    beta[t + 1, ])
      
      xi[t,,] = xi[t,,] - matrixStats::logSumExp(xi[t,,])
    }
    
    # M step
    # emission
    denominators = matrixStats::colLogSumExps(gamma)
    exp_gamma_denominators = exp(t(t(gamma) - denominators))
    for (i in seq_len(hmm_param$num_features)) {
      emission_paras[["mu"]][, i] = 
        exp(matrixStats::colLogSumExps(gamma + log(observation[, i])) - denominators)
    }
    
    sd_tmp = emission_paras[["sd"]]
    for (i in seq_len(hmm_param$num_features)) {
      for (j in states) {
        sd_tmp[j, i] = 
          exp((matrixStats::logSumExp(gamma[, j] +
                                        log((observation[, i] -
                                               emission_paras[["mu"]][j, i])^2)
                                      ) - denominators[j]) / 2)
      }
    }
    # limit min sd
    sd_tmp = apply(sd_tmp, 2, 
                   function(x) {
                     x[x < 0.1 | is.na(x)] = median(x, na.rm = TRUE)
                     x })
    emission_paras[["sd"]] = sd_tmp
    
    # transition
    denominators = matrixStats::colLogSumExps(gamma[-.len, ])
    
    for (i in states) { 
      numerators = matrixStats::colLogSumExps(xi[-.len, i, ])
      
      log_transition_probs[i, ] = numerators - denominators[i]
      
    }
    
    # update parameters
    hmm_param$initstate_probs = exp(gamma[1, ]) + 1e-20
    hmm_param$log_emission_probs = get_emission_probs(observation, emission_paras)
    hmm_param$emission_paras = emission_paras
    hmm_param$log_transition_probs = log_transition_probs
    
    hmm_param
  }
  
  
  hmm_iterator <- function(observation, hmm_param, 
                           max_iter = 20, show_plot = TRUE,
                           tolerence = 3) {
    
    hmm_param_update = hmm_param
    
    alpha_end_probs = c()
    
    for (i in seq_len(max_iter)) {
      
      # likelihood
      alpha_end_probs = 
        c(alpha_end_probs, 
          matrixStats::logSumExp(tail(forward(observation, hmm_param_update), 1)) )
      message(paste0("Iteration ", i, ": log-likelihood ", alpha_end_probs[i]))
      
      # stop if increment is small compared with the initial states
      if (i > 1) {
        if (diff(alpha_end_probs[c(i-1, i)]) / diff(alpha_end_probs[c(1, i)]) < 1e-5) {
          tolerence = tolerence - 1
          if (tolerence == 0) break
        }
      }
      
      hmm_param_update = baum_welch(observation, hmm_param_update)
      
      if (show_plot) {
        # png(paste0("HMM_emission_state_mean_iteration_", i, ".png"), width = 800, height = 500)
        mu = hmm_param_update$emission_paras$mu
        mu = t(mu) / colMeans(mu, na.rm = TRUE)
        image(mu[, rev(seq_len(ncol(mu)))], main = paste("Iteration", i),
              axes = FALSE, ylab = "States")
        axis(side = 1, at = seq(0, 1, length.out = ncol(observation)), labels = NA)
        text(x = seq(0, 1, length.out = ncol(observation)),
             y = par("usr")[3] - 0.1,
             labels = colnames(observation),
             xpd = NA,
             adj = 0.8,
             srt = 35)
        axis(side = 2, at = seq(0, 1, length.out = length(hmm_param$states)), 
             labels = rev(hmm_param$states), las = 2)
        # dev.off()
      }
      
    }
    
    hmm_param_update$alpha_end_probs = alpha_end_probs
    hmm_param_update
  }
  
  #  ----------- run -----------
  # adjust feature negative scales if exist
  feature_mins <- colMins(observation)
  observation <- t(t(observation) - feature_mins)
  
  # initiate
  hmm_param <- initializer(num_states = num_states, observation, init_method = "kmeans")
  
  if (F) {  # minor steps
    alpha = forward(observation, hmm_param)
    head(alpha)
    tail(alpha)
    
    beta = backward(observation, hmm_param)
    head(beta)
    tail(beta)
    
    hmm_param <- baum_welch(observation, hmm_param)
    hmm_param$emission_paras$mu
    hmm_param$emission_paras$sd
    head(t(hmm_param$log_emission_probs))
    tail(t(hmm_param$log_emission_probs))
  }
    
  # main steps
  hmm_param <- hmm_iterator(observation, hmm_param, max_iter = 40, tolerence = 5)
  hmm_param$emission_paras$mu <- t(t(hmm_param$emission_paras$mu) + feature_mins)
  # print(hmm_param$emission_paras$mu)
  
  # make prediction
  viterbi_path <- viterbi(observation, hmm_param)[[2]]
  table(viterbi_path)
  
  gc()
  
  # output
  message("Done.")
  list("hmm_param" = hmm_param,
       "viterbi_path" = viterbi_path)
  
} # end of the wrapper 
