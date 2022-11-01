# Rui Shao, 2022 Oct
#
# -----------------------------------------------------------------------------
#
# Genome segmentation with continuous inputs
# Hidden Markov model, Baum-Welch algorithm


# --------------------------------- Steps -------------------------------------
#
# 0. hmm_initializer()        # random log-normal centroids
#
# 1. get_emission_probs()
#
# 2. baum_welch() 
#   |
#   |---- forward() and backward() to get alpha and beta
#   |
#   `---- gamma and xi to EM update transition_prob and emission_prob
#
# 3. viterbi(), calculate the optimal state path with viterbi method
#
# -----------------------------------------------------------------------------


# ------------------------------ Functions ------------------------------------
#
# 1. get_emission_probs()
#
# 2. hmm_initializer()
#
# 3. forward() 
#
# 4. backward()
#
# 5. baum_welch()
#
# 6. get_viterbi()
#
# 7. hmm_iterator()
#
# 8. show_state_mu()
#
# 9. run_hmm()
#
# -----------------------------------------------------------------------------


# ------------------------------ Approaches -----------------------------------
#
# i. concatenated
#   input: list(A, B, C)
# 
# ii. averaged
#   input: list((A + B + C) / 3)
#
# -----------------------------------------------------------------------------

# Cpp functions
Rcpp::cppFunction(
  "NumericMatrix get_alpha_Cpp(NumericMatrix alpha, 
                                 NumericVector log_initstate_probs,
                                 NumericMatrix log_transition_probs, 
                                 NumericMatrix log_emission_probs
                                 ) {
                int n = alpha.nrow(), k = alpha.ncol();
                
                for (int j = 0; j < k; ++j) {
                  alpha(0, j) = log_initstate_probs(j) + log_emission_probs(j, 0);
                }
                
                for (int i = 1; i < n; ++i) 
                {
                  for (int j = 0; j < k; ++j) 
                  {
                    double tmp = alpha(i-1, 0) + 
                      log_transition_probs(0, j) + 
                      log_emission_probs(j, i);
                    
                    for (int x = 1; x < k; ++x) 
                    {
                      double tmp_next = alpha(i-1, x) + 
                      log_transition_probs(x, j) + 
                      log_emission_probs(j, i);
                      
                      if (tmp_next > tmp) 
                      {
                        tmp = tmp_next + log(1 + exp(tmp - tmp_next));
                      } else 
                      {
                        tmp = tmp + log(1 + exp(tmp_next - tmp));
                      }
                    }
                    
                    alpha(i, j) = tmp;
                  }
                }
                return alpha;
              }")

Rcpp::cppFunction(
  "NumericMatrix get_beta_Cpp(NumericMatrix beta, 
                              NumericMatrix log_transition_probs, 
                              NumericMatrix log_emission_probs 
                              ) {
                int n = beta.nrow(), k = beta.ncol();
  
                for (int i = n - 2; i >= 0; i--) 
                {
                  for (int j = 0; j < k; ++j) 
                  {
                    double tmp = beta(i + 1, 0) + 
                      log_transition_probs(j, 0) + 
                      log_emission_probs(0, i + 1);
                    
                    for (int x = 1; x < k; ++x) 
                    {
                      double tmp_next = beta(i + 1, x) + 
                      log_transition_probs(j, x) + 
                      log_emission_probs(x, i + 1);
                      
                      if (tmp_next > tmp) 
                      {
                        tmp = tmp_next + log(1 + exp(tmp - tmp_next));
                      } else 
                      {
                        tmp = tmp + log(1 + exp(tmp_next - tmp));
                      }
                    }
                    beta(i, j) = tmp;
                  }
                }
                return beta;
              }")

Rcpp::cppFunction(
  'IntegerVector get_viterbi_Cpp(List hmm_param,
                                 List viterbi_paths,
                                 int num_obs
                                 ) {
                int n = num_obs, k = hmm_param["num_states"];
                
                NumericMatrix log_emission_probs = hmm_param["log_emission_probs"];
                NumericMatrix log_transition_probs = hmm_param["log_transition_probs"];
                NumericVector initstate_probs = hmm_param["initstate_probs"];
                
                NumericVector viterbi_probs = log(initstate_probs) + log_emission_probs(_, 0);
                
                for (int i = 1; i < n; ++i) 
                {
                  // List tmp_viterbi_paths = viterbi_paths;
                  NumericVector viterbi_tmp (k);
                  
                  for (int x = 0; x < k; ++x) 
                  {
                    NumericVector intermediate_probs = viterbi_probs + log_transition_probs(_, x);
                    
                    // find the max state
                    int max_state = 0;
                    double max_prob = intermediate_probs(0);
                    for (int y = 0; y < k; ++y) 
                    {
                      if (max_prob < intermediate_probs(y)) 
                      {
                        max_state = y;
                        max_prob = intermediate_probs(y);
                      }
                    }
                    
                    viterbi_tmp(x) = log_emission_probs(x, i) + max_prob; 
                    
                    IntegerVector tmp_max_path = viterbi_paths[max_state]; 
                    tmp_max_path = clone(tmp_max_path);
                    tmp_max_path[i] = x + 1;
                    viterbi_paths[x] = tmp_max_path; 
                  }
                  viterbi_probs = viterbi_tmp; 
                }
                
                // find which max state
                int max_path = 0;
                double max_prob = viterbi_probs(0);
                for (int y = 0; y < k; ++y) 
                {
                  if (max_prob < viterbi_probs(y)) 
                  {
                    max_path = y;
                    max_prob = viterbi_probs(y);
                  }
                }
                return viterbi_paths[max_path];
              }')

Rcpp::cppFunction(
  "NumericVector get_xi_Cpp(NumericVector xi, 
                            NumericMatrix alpha, 
                            NumericMatrix beta, 
                            NumericMatrix log_transition_probs, 
                            NumericMatrix log_emission_probs 
                            ) {
                int n = alpha.nrow(), k = alpha.ncol();
                
                for (int i = 0; i < n - 1; ++i) 
                {
                  double prob_sum = R_NegInf;
                  for (int j = 0; j < k; ++j) //state from
                  {
                    for (int x = 0; x < k; ++x) //state to
                    {
                      double prob = alpha(i, j) + 
                      log_transition_probs(j, x) + 
                      log_emission_probs(x, i + 1) +
                      beta(i + 1, x);
                      
                      xi(j + x * k + i * k * k) = prob;
                      
                      if (j == 0 && x == 0) 
                      {
                        prob_sum = prob;
                      } 
                      else if (prob > prob_sum) 
                      {
                        prob_sum = prob + log(1 + exp(prob_sum - prob));
                      } 
                      else 
                      {
                        prob_sum = prob_sum + log(1 + exp(prob - prob_sum));
                      }
                    }
                  }
                  
                  for (int j = 0; j < k; ++j) 
                  {
                    for (int x = 0; x < k; ++x) 
                    {
                      xi(j + x * k + i * k * k) -= prob_sum;
                    }
                  }
                }
                return xi;
              }")

Rcpp::cppFunction(
  'IntegerVector get_viterbi_Cpp(List hmm_param,
                                 List viterbi_paths,
                                 int num_obs
                                 ) {
                int n = num_obs, k = hmm_param["num_states"];
                
                NumericMatrix log_emission_probs = hmm_param["log_emission_probs"];
                NumericMatrix log_transition_probs = hmm_param["log_transition_probs"];
                NumericVector initstate_probs = hmm_param["initstate_probs"];
                
                NumericVector viterbi_probs = log(initstate_probs) + log_emission_probs(_, 0);
                
                for (int i = 1; i < n; ++i) 
                {
                  // List tmp_viterbi_paths = viterbi_paths;
                  NumericVector viterbi_tmp (k);
                  
                  for (int x = 0; x < k; ++x) 
                  {
                    NumericVector intermediate_probs = viterbi_probs + log_transition_probs(_, x);
                    
                    // find the max state
                    int max_state = 0;
                    double max_prob = intermediate_probs(0);
                    for (int y = 0; y < k; ++y) 
                    {
                      if (max_prob < intermediate_probs(y)) 
                      {
                        max_state = y;
                        max_prob = intermediate_probs(y);
                      }
                    }
                    
                    viterbi_tmp(x) = log_emission_probs(x, i) + max_prob; 
                    
                    IntegerVector tmp_max_path = viterbi_paths[max_state]; 
                    tmp_max_path = clone(tmp_max_path);
                    tmp_max_path[i] = x + 1;
                    viterbi_paths[x] = tmp_max_path; 
                  }
                  viterbi_probs = viterbi_tmp; 
                }
                
                // find which max state
                int max_path = 0;
                double max_prob = viterbi_probs(0);
                for (int y = 0; y < k; ++y) 
                {
                  if (max_prob < viterbi_probs(y)) 
                  {
                    max_path = y;
                    max_prob = viterbi_probs(y);
                  }
                }
                return viterbi_paths[max_path];
              }')

#' Baum-welch HMM with log-normal probability estimation
#' @param observations A list of matrices with features in normal distribution
#' @param num_states  Number of hidden states
#' @param init_method Initiation HMM states with kmeans, mclust, or quantile methods.
#' 
hmm_initializer <- function(observations = NULL, 
                            num_states  = 10, 
                            init_method = c("kmeans", "mclust", "quantile"),
                            seed = 1) {
  # Args:
  # num_states:  number of hidden states
  # observation: input observation matrix
  
  stopifnot(!is.null(observations) & length(observations) > 0)
  set.seed(seed)
  
  num_features = ncol(observations[[1]])
  num_obsers = unlist(lapply(observations, nrow))
  states = seq_len(num_states)
  
  # extract 10% data for initialization
  num_obsers_small = num_obsers %/% 10
  num_obsers_small_start = cumsum(num_obsers_small) - num_obsers_small + 1
  observation = matrix(0, nrow = sum(num_obsers_small), ncol = num_features)
  for (i in seq_along(observations)) {
    observation[seq(num_obsers_small_start[i], cumsum(num_obsers_small)[i]), ] =
      observations[[i]][sample(num_obsers[i], num_obsers_small[i]), ]
  }
  
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
  if (!is.null(colnames(observations[[1]]))) 
    colnames(initi_state_mean) = colnames(observations[[1]])
  
  # set initial state sd at total sd level
  initi_state_sd = apply(observation, 2, 
                         function(x) 
                           rep(sd(x) / num_states, num_states) 
  ) 
  rownames(initi_state_sd) = rownames(initi_state_mean)
  colnames(initi_state_sd) = colnames(initi_state_mean)
  
  # states X observs
  emission_paras = list("mu" = initi_state_mean,  # mean of each state
                        "sd" = initi_state_sd)  # standard deviation of each state
  # log_emission_probs = get_emission_probs(observation, emission_paras) # log probs
  
  transition_probs = matrix(1, nrow = num_states, ncol = num_states) / num_states
  dimnames(transition_probs) = list(states, states)
  
  # create a parameter list 
  hmm_param = list()
  hmm_param$states = states
  hmm_param$num_states = num_states
  hmm_param$num_features = num_features
  
  hmm_param$initi_state_mean = initi_state_mean
  hmm_param$initi_state_sd = initi_state_sd
  
  hmm_param$initstate_probs = rep(1, num_states) / num_states
  hmm_param$log_transition_probs = log(transition_probs)
  
  hmm_param$emission_paras = emission_paras
  # hmm_param$log_emission_probs = log_emission_probs
  hmm_param
}


#' Multi-normal distribution for emission probabilities
#' 
get_emission_probs <- function(observation, emission_paras) {
  # output: (state, position) emission matrix in log
  num_states = nrow(emission_paras[["mu"]])
  log_emission_probs <- matrix(0, nrow = num_states, ncol = nrow(observation))
  
  for (i in seq_len(num_states)) { # joint probability
    log_emission_probs[i, ] = 
      colSums(log(dnorm(t(observation), 
                        mean = emission_paras[["mu"]][i, ], 
                        sd = emission_paras[["sd"]][i, ]) + 1e-4))
  }
  log_emission_probs[log_emission_probs > 0] = 0
  rownames(log_emission_probs) <- rownames(emission_paras[["mu"]])
  log_emission_probs
}



#' HMM parameters trainning with Baum-Welch
#' @param observation
#' @param hmm_param
#' 
# forward 
forward <- function(observation, hmm_param, is_last = FALSE) {
  # Args:
  # observation: position X feature matrix
  # hmm_param:   a parameter list
  
  alpha = matrix(0, ncol = hmm_param$num_states, nrow = nrow(observation))
  colnames(alpha) = hmm_param$states
  log_emission_probs = get_emission_probs(observation, hmm_param$emission_paras)
  
  # output in log
  alpha = get_alpha_Cpp(alpha = alpha, 
                        log_initstate_probs = log(hmm_param$initstate_probs), 
                        log_transition_probs = hmm_param$log_transition_probs, 
                        log_emission_probs = log_emission_probs)
  
  if (is_last)  {
    return(tail(alpha, 1))
  } else {
    return(alpha)
  }
}


# backward
backward <- function(observation, hmm_param) {
  # Args:
  # observation: position X feature matrix
  # hmm_param:   a parameter list
  
  beta = matrix(0, ncol = hmm_param$num_states, nrow = nrow(observation))
  colnames(beta) = hmm_param$states
  log_emission_probs = get_emission_probs(observation, hmm_param$emission_paras)
  
  # output in log
  get_beta_Cpp(beta = beta, 
               log_transition_probs = hmm_param$log_transition_probs, 
               log_emission_probs = log_emission_probs)
}


# update parameters
baum_welch <- function(observation, hmm_param) {
  # Args:
  # observation: position X feature matrix
  # hmm_param:   a parameter list
  
  if (nrow(observation) == 0) return(NULL)
  
  # start 
  alpha = forward(observation = observation, hmm_param = hmm_param, is_last = FALSE)
  beta = backward(observation = observation, hmm_param = hmm_param)  
  log_emission_probs = get_emission_probs(observation, hmm_param$emission_paras)
  
  xi = array(0, 
             dim = c(hmm_param$num_states,
                     hmm_param$num_states, 
                     nrow(observation) - 1), 
             dimnames = list(hmm_param$states, 
                             hmm_param$states, 
                             seq_len(nrow(observation) - 1)))
  
  gamma = matrix(0, 
                 ncol = hmm_param$num_states, 
                 nrow = nrow(observation))
  colnames(gamma) = hmm_param$states
  
  # E step
  xi = get_xi_Cpp(xi = xi,
                  alpha = alpha,
                  beta = beta, 
                  log_transition_probs = hmm_param$log_transition_probs,
                  log_emission_probs = log_emission_probs)
  
  gamma = alpha + beta
  gamma = gamma - matrixStats::rowLogSumExps(gamma)
  
  # M step
  # update initial parameters
  hmm_param$initstate_probs = exp(gamma[1, ])
  
  # update emission
  denominators = matrixStats::colLogSumExps(gamma)
  exp_gamma_denominators = exp(t(t(gamma) - denominators))
  
  # update mu 
  for (i in seq_len(hmm_param$num_features)) {
    hmm_param$emission_paras[["mu"]][, i] = 
      exp(matrixStats::colLogSumExps(gamma + log(observation[, i])) - denominators)
  }
  
  # update sd
  sd_tmp = hmm_param$emission_paras[["sd"]]
  for (i in seq_len(hmm_param$num_features)) {
    for (j in hmm_param$states) {
      sd_tmp[j, i] = 
        exp((matrixStats::logSumExp(gamma[, j] +
                                      log((observation[, i] -
                                             hmm_param$emission_paras[["mu"]][j, i])^2)
        ) - denominators[j]) / 2)
    }
    # limit min sd
    # sd_tmp_i = sd_tmp[, i]
    # idx_i = is.na(sd_tmp_i) | sd_tmp_i < (hmm_param$initi_state_sd[, i] / hmm_param$num_states / 5)
    # sd_tmp_i[idx_i] = hmm_param$initi_state_sd[idx_i, i] / sqrt(hmm_param$num_states)
    # sd_tmp[, i] = sd_tmp_i
  }
  hmm_param$emission_paras[["sd"]] = sd_tmp
  
  # update transition
  denominators = matrixStats::colLogSumExps(gamma[-nrow(observation), ])
  
  for (i in hmm_param$states) { 
    numerators = matrixStats::rowLogSumExps(xi[i, , -nrow(observation)])
    hmm_param$log_transition_probs[i, ] = numerators - denominators[i]
  }
  
  #
  # hmm_param$log_transition_probs[hmm_param$log_transition_probs < -25] = -25
  # hmm_param$log_transition_probs[hmm_param$log_transition_probs > 0] = 0
  
  # hmm_param$log_emission_probs = get_emission_probs(observation, 
  #                                                   hmm_param$emission_paras)
  hmm_param
}


average_parameters <- function(hmm_param_list) {
  hmm_param = hmm_param_list[[1]]
  for (i in seq_along(hmm_param_list)[-1]) {
    hmm_param$initstate_probs = rbind(hmm_param$initstate_probs, 
                                      hmm_param_list[[i]]$initstate_probs)
    hmm_param$log_transition_probs = abind::abind(hmm_param$log_transition_probs, 
                                                  hmm_param_list[[i]]$log_transition_probs, along = 3)
    hmm_param$emission_paras$mu = abind::abind(hmm_param$emission_paras$mu, 
                                                  hmm_param_list[[i]]$emission_paras$mu, along = 3)
    hmm_param$emission_paras$sd = abind::abind(hmm_param$emission_paras$sd, 
                                        hmm_param_list[[i]]$emission_paras$sd, along = 3)
  }
  
  hmm_param$log_transition_probs = apply(hmm_param$log_transition_probs, 2, 
                                         function(x) {
                                           matrixStats::rowLogSumExps(x) - log(ncol(x))
                                         })
  hmm_param$initstate_probs = colMeans(hmm_param$initstate_probs)
  hmm_param$emission_paras$mu = apply(hmm_param$emission_paras$mu, 2, 
                                      function(x) rowMeans(x))
  hmm_param$emission_paras$sd = apply(hmm_param$emission_paras$sd, 2, 
                                      function(x) rowMeans(x))
  
  hmm_param
}


hmm_iterator <- function(observation_list,
                         hmm_param,
                         max_iter = 20, 
                         tolerence = 5,
                         show_plot = FALSE, 
                         seed = 1) {
  require(foreach)
  require(doParallel)
  registerDoParallel(cores = 4)
  set.seed(seed)
  n_obs = length(observation_list)
  cat("Process", n_obs, ifelse(n_obs > 1, "sequences.\n", "sequence.\n"))
  
  alpha_probs = c()
  
  for (iter in seq_len(max_iter)) {
    
    # update emission and transition parameters
    hmm_param_list = list()
    
    # each chromosome as a small batch will cause transition probs unstable
    hmm_param_list = foreach(obs = observation_list) %dopar% {
      baum_welch(obs, hmm_param)
    }
    
    # for (i in seq_len(n_obs)) {
    #   cat(paste("\r", (i * 100) %/% n_obs),
    #       '% |',
    #       rep('=', i * 50 / n_obs),
    #       ifelse(i == n_obs, "|\n", ">"),
    #       sep = '')
    #   hmm_param_list = c(hmm_param_list,
    #                      list(baum_welch(observation_list[[n_obs_rand[i]]], hmm_param)))
    # }
    hmm_param = average_parameters(hmm_param_list)
    
    # compute new alpha
    alpha_last_mat = foreach(obs = observation_list, .combine = rbind) %dopar% {
      forward(obs, hmm_param, is_last = TRUE)
    }
    
    # forward likelihood
    alpha_probs = c(alpha_probs, matrixStats::logSumExp(colSums(alpha_last_mat)))
    message(paste0("Iteration ", iter, ": log-likelihood ", alpha_probs[iter]))
    
    # stop if increment is small compared with the initial states
    if (iter > 1) {
      if (any(is.na(alpha_probs))) break
      if (diff(alpha_probs[c(iter - 1, iter)]) / diff(alpha_probs[c(1, iter)]) < 1e-5) {
        tolerence = tolerence - 1
        if (tolerence == 0) {
          message("Parameter training finished.")
          break
        }
      }
    }
    
    if (show_plot) show_state_mu(hmm_param, paste("Iteration", iter))
  }
  
  hmm_param
}



show_state_mu <- function(hmm_param, num_features, col_names, .title = "State averages") {
  # hmm_param: plot emission mu
  # iter: number of iterations
  mu = hmm_param$emission_paras$mu
  mu = t(mu) / colMeans(mu, na.rm = TRUE)
  image(mu[, rev(seq_len(ncol(mu)))], main = .title,
        axes = FALSE, ylab = "States")
  axis(side = 1, at = seq(0, 1, length.out = num_features), labels = NA)
  text(x = seq(0, 1, length.out = num_features),
       y = par("usr")[3] - 0.1,
       labels = col_names,
       xpd = NA,
       adj = 0.8,
       srt = 35)
  axis(side = 2, at = seq(0, 1, length.out = length(hmm_param$states)),
       labels = rev(hmm_param$states), las = 2)
}



#' Compute viterbi path after training HMM parameters  
#' @param observation
#' @param hmm_param
#' 
get_viterbi <- function(observation_list, hmm_param) {
  
  require(foreach)
  require(doParallel)
  registerDoParallel(cores = 4)
  
  message("Get viterbi paths.")
  # get viterbi paths for each chromosome
  
  viterbi_path_list = foreach(obs = observation_list) %dopar% {
    hmm_param$log_emission_probs = get_emission_probs(observation = obs,
                                                      emission_paras = hmm_param$emission_paras)
    init_viterbi_paths = lapply(seq_len(hmm_param$num_states), 
                                function(x) as.integer(c(x, numeric(nrow(obs) - 1))))
    get_viterbi_Cpp(hmm_param,
                    init_viterbi_paths, 
                    nrow(obs))
  }
  
  if (!is.null(names(observation_list))) names(viterbi_path_list) = names(observation_list)
  viterbi_path_list
}


run_hmm <- function(observation_list,
                    num_states = 10, 
                    init_method = "kmeans",
                    seed = 1) {
  
  # --------------------- run ------------------------
  # # adjust feature negative scales if exist
  # feature_mins <- matrixStats::colMins(observation)
  # observation <- t(t(observation) - feature_mins)
  
  # initiate
  
  hmm_param <- hmm_initializer(observations = observation_list,
                               num_states = num_states, 
                               init_method = init_method,
                               seed = seed)
  
  if (TRUE) {  # initialization steps
    alpha = forward(observation_list[[1]], hmm_param)
    message("Forward ok.")
    
    beta = backward(observation_list[[1]], hmm_param)
    message("Backward ok.")
    
    hmm_param <- baum_welch(observation_list[[1]], hmm_param)
    message("Baum-welch ok.\n")
  }
  
  # main steps
  hmm_param <- hmm_iterator(observation_list, hmm_param, max_iter = 50, seed = seed)
  # hmm_param$emission_paras$mu <- t(t(hmm_param$emission_paras$mu) + feature_mins)
  show_state_mu(hmm_param,
                ncol(observation_list[[1]]), 
                colnames(observation_list[[1]]))
  
  # make prediction
  viterbi_paths <- get_viterbi(observation_list, hmm_param)
  print(sapply(viterbi_paths, table))
  
  # output
  list("hmm_param" = hmm_param,
       "viterbi_paths" = viterbi_paths)
  
} # end of the wrapper 
