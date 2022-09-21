# Rui Shao, 2022 Sep
# Baum-Welch algorithm for discrete Hidden Markov models

# -----------------------------------------------------------------------------

# DIY for genomic region segmentation

# -----------------------------------------------------------------------------

run_hmm <- function(observation, num_states = 10) {
  # Args:
  # num_states: number of hidden states
  # observation: input observation matrix, log-transformed
  
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
    
    # output: (state, position) matrix
    
    num_states = nrow(emission_paras[["mu"]])
    emission_probs <- matrix(1, nrow = num_states, ncol = nrow(observation))
    
    for (i in seq_len(num_states)) { # joint probability
      emission_probs[i, ] = 
        matrixStats::colProds(dnorm(t(observation), 
                                    mean = emission_paras[["mu"]][i, ], 
                                    sd = emission_paras[["sd"]][i, ]))
    }
    rownames(emission_probs) <- rownames(emission_paras[["mu"]])
    emission_probs
  }
  
  
  initializer <- function(num_states = 10,  observation) {
    # Args:
    # num_states: number of hidden states
    # observation: input observation matrix, log-transformed
    
    num_features = ncol(observation)
    states = seq_len(num_states)
    
    # sample initial state means
    initi_state_mean = apply(observation, 2, 
                             function(x) 
                               sample(quantile(x, seq(0, 1, length.out = num_states)))
    )
    rownames(initi_state_mean) = states
    
    initi_state_mu = apply(observation, 2, 
                           function(x) 
                             rep(sd(x), num_states) 
    ) 
    rownames(initi_state_mu) = states
    
    # states X observs
    emission_paras = list("mu" = initi_state_mean,  # mean of each state
                          "sd" = initi_state_mu)  # standard deviation of each state
    emission_probs = get_emission_probs(observation, hmm_param$emission_paras)
    
    transition_probs = matrix(1, nrow = num_states, ncol = num_states) / num_states
    dimnames(transition_probs) = list(states, states)
    
    # create a parameter list 
    hmm_param = list()
    hmm_param$states = states
    hmm_param$num_features = num_features
    hmm_param$initstate_probs = rep(1, num_states) / num_states
    hmm_param$transition_probs = transition_probs
    hmm_param$emission_paras = emission_paras
    hmm_param$emission_probs = emission_probs
    hmm_param
  }
  
  
  forward <- function(observation, hmm_param) {
    # Args:
    # observation: position X feature matrix
    # hmm_param: a parameter list
    
    states = hmm_param$states
    log_initstate_probs  = log(hmm_param$initstate_probs)
    log_transition_probs = log(hmm_param$transition_probs)
    log_emission_probs = log(hmm_param$emission_probs)
    
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
    log_emission_probs   = log(hmm_param$emission_probs)
    log_transition_probs = log(hmm_param$transition_probs)
    
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
    
    states =  hmm_param$states
    log_initstate_probs  = log(hmm_param$initstate_probs)
    log_emission_probs   = log(hmm_param$emission_probs)
    log_transition_probs = log(hmm_param$transition_probs)
    
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
    
    max_path = factor(state_path[[which.max(viterbi)]],
                      levels = states[order(rowMeans(hmm_param$emission_paras$mu),
                                            decreasing = TRUE)])
    
    
    list(max(viterbi), max_path)
  }
  
  
  baum_welch <- function(observation, hmm_param) {
    # Args:
    # observation: position X feature matrix
    # hmm_param: a parameter list
    
    states =  hmm_param$states
    emission_paras = hmm_param$emission_paras
    log_emission_probs   = log(hmm_param$emission_probs)
    transition_probs = hmm_param$transition_probs
    log_transition_probs = log(transition_probs)
    
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
    
    for (i in seq_len(hmm_param$num_features)) {
      emission_paras[["mu"]][, i] = 
        exp(matrixStats::colLogSumExps(gamma + log(observation[, i])) - denominators)
    }
    
    for (i in seq_len(hmm_param$num_features)) {
      for (j in states) {
        emission_paras[["sd"]][j, i] = 
          exp((matrixStats::logSumExp(gamma[, j] + 
                                        log((observation[, i] - 
                                               emission_paras[["mu"]][j, i])^2) 
                                      ) - denominators[j]) / 2)
      }
    }
    emission_paras[["sd"]][emission_paras[["sd"]] == 0] = 
      quantile(emission_paras[["sd"]], 0.1) # limit zero
    
    # transition
    denominators = matrixStats::colLogSumExps(gamma[-.len, ])
    
    for (i in states) { 
      numerators = matrixStats::colLogSumExps(xi[-.len, i, ])
      
      transition_probs[i, ] = exp(numerators - denominators[i])
      
    }
    transition_probs[is.na(transition_probs)] = 1e-200
    # print(transition_probs)
    
    
    # update parameters
    hmm_param$initstate_probs = exp(gamma[1, ]) + 1e-200
    hmm_param$emission_probs = get_emission_probs(observation, emission_paras)
    hmm_param$emission_paras = emission_paras
    hmm_param$transition_probs = transition_probs
    
    hmm_param
  }
  
  
  hmm_iterator <- function(observation, hmm_param, max_iter = 20, show_plot = TRUE) {
    
    hmm_param_update = hmm_param
    
    alpha_end_probs = c()
    
    for (i in seq_len(20)) {
      
      alpha_end_probs = 
        c(alpha_end_probs, 
          matrixStats::logSumExp(tail(forward(observation, hmm_param_update), 1)) )
      
      # stop if increment is smaller than 0.001 %
      if (i > 1) {
        if (diff(alpha_end_probs[c(i-1, i)]) / diff(alpha_end_probs[c(1, i)]) < 1e-6) break
      }
      
      hmm_param_update = baum_welch(observation, hmm_param_update)
      
      if (show_plot)
        image(t(hmm_param_update$emission_paras$mu), main = paste("Iteration", i))
    }
    
    hmm_param_update$alpha_end_probs = alpha_end_probs
    hmm_param_update
  }
  
  # run -----------
  
  hmm_param <- initializer(num_states = 10, observation)
  
  hmm_param <- hmm_iterator(observation, hmm_param)
  
  viterbi_path <- viterbi(observation, hmm_param)[[2]]

  list("hmm_param" = hmm_param,
       "viterbi_path" = viterbi_path)
  
} # end of wrapper 
