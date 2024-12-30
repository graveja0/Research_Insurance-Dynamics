#' Extract transition rates from MSM fit
#' @param model Fitted msm model
#' @param ages Vector of ages at which to extract rates
get_transition_rates <- function(model, ages) {
    # Get rates at each age
    rates_by_age <- map_dfr(ages, function(age) {
        Q <- get_transition_intensities(model, age)
        tibble(
            age = age,
            from_state = rep(rownames(Q), each = ncol(Q)),
            to_state = rep(colnames(Q), times = nrow(Q)),
            rate = as.vector(Q)
        ) %>%
            filter(from_state != to_state)  # Remove diagonal elements
    })
    return(rates_by_age)
}

#' Create flexible parametric functions for transition rates
#' @param knots Knot locations for splines
create_rate_functions <- function(knots) {
    # Create basis functions
    basis <- ns(knots, df = length(knots)-2)
    
    # Return function that generates rates given parameters
    function(age, params, from_state, to_state) {
        # Get spline values for this age
        spline_vals <- predict(basis, age)
        # Calculate log rate (ensure positivity)
        log_rate <- sum(params * spline_vals)
        return(exp(log_rate))
    }
}

#' Simulate system forward using given rate functions
#' @param rate_functions List of rate functions
#' @param params Parameters for rate functions
#' @param initial_dist Initial distribution
simulate_system <- function(rate_functions, params, initial_dist, ages) {
    n_states <- length(initial_dist)
    results <- matrix(0, nrow = length(ages), ncol = n_states)
    current_dist <- initial_dist
    
    # Print diagnostics
    cat("Number of states:", n_states, "\n")
    cat("Rate function names:", names(rate_functions), "\n")
    cat("Parameter names:", names(params), "\n")
    
    for(i in seq_along(ages)) {
        age <- ages[i]
        results[i,] <- current_dist
        
        # Build Q matrix for this age
        Q <- matrix(0, n_states, n_states)
        for(from in 1:n_states) {
            for(to in 1:n_states) {
                if(from != to) {
                    transition_name <- paste(from, to, sep="-")
                    if(transition_name %in% names(rate_functions)) {
                        rate_fn <- rate_functions[[transition_name]]
                        Q[from,to] <- rate_fn(
                            age = age,
                            params = params[[transition_name]],
                            from_state = from,
                            to_state = to
                        )
                    }
                }
            }
            Q[from,from] <- -sum(Q[from,])
        }
        
        # Update distribution
        P <- expm(Q)
        current_dist <- current_dist %*% P
    }
    return(results)
}
#' Calibrate using Incremental Mixture Importance Sampling
#' @param model Initial MSM model
#' @param targets Target occupancy data
#' @param n_particles Number of particles for IMIS
calibrate_transitions <- function(model, targets, n_particles = 1000) {
    require(IMIS)
    
    # Setup transitions and parameters
    states <- attr(model, "states")
    transitions <- expand.grid(
        from = seq_along(states),
        to = seq_along(states)
    ) %>%
        filter(from != to) %>%
        mutate(
            from_state = states[from],
            to_state = states[to],
            name = paste(from, to, sep="-")
        )
    
    # Setup parameters
    knots <- attr(model, "knots")
    rate_functions <- map(transitions$name, ~create_rate_functions(knots))
    names(rate_functions) <- transitions$name
    n_params_per_transition <- ncol(ns(knots, df = length(knots)-2))
    n_transitions <- nrow(transitions)
    total_params <- n_params_per_transition * n_transitions
    
    cat("Testing basic functionality:\n")
    cat("Number of states:", length(states), "\n")
    cat("Number of transitions:", n_transitions, "\n")
    cat("Parameters per transition:", n_params_per_transition, "\n")
    cat("Total parameters:", total_params, "\n\n")
    
    # Test likelihood calculation with different parameter values
    test_likelihood <- function(scale) {
        test_params <- rnorm(total_params, mean = 0, sd = scale)
        params <- split(test_params, rep(1:n_transitions, each = n_params_per_transition))
        names(params) <- transitions$name
        
        init_dist <- targets %>% 
            filter(agep == min(agep)) %>% 
            pull(pct)
        
        cat("Testing with scale", scale, ":\n")
        cat("Initial distribution:", round(init_dist, 3), "\n")
        
        predicted <- simulate_system(rate_functions, params, init_dist, unique(targets$agep))
        cat("Predicted range:", range(predicted), "\n")
        
        ll <- sum(dnorm(
            x = as.vector(predicted),
            mean = targets$pct,
            sd = targets$pct_se,
            log = TRUE
        ))
        cat("Log-likelihood:", ll, "\n\n")
        
        return(ll)
    }
    
    # Test with different scales
    scales <- c(0.001, 0.01, 0.1)
    test_results <- sapply(scales, test_likelihood)
    
    cat("Test results for different scales:\n")
    print(test_results)
    
    # If we get here without errors, proceed with actual calibration
    # Use the scale that gave the best likelihood
    best_scale <- scales[which.max(test_results)]
    cat("\nProceeding with scale:", best_scale, "\n")
    
    # Define IMIS functions using the best scale
    sample.prior <- function(n) {
        matrix(rnorm(n * total_params, mean = 0, sd = best_scale),
               nrow = n, 
               ncol = total_params)
    }
    
    prior <- function(theta) {
        if(is.vector(theta)) theta <- matrix(theta, nrow = 1)
        apply(theta, 1, function(x) {
            sum(dnorm(x, mean = 0, sd = best_scale, log = TRUE))
        })
    }
    
    likelihood <- function(theta) {
        if(is.vector(theta)) theta <- matrix(theta, nrow = 1)
        apply(theta, 1, function(x) {
            params <- split(x, rep(1:n_transitions, each = n_params_per_transition))
            names(params) <- transitions$name
            predicted <- simulate_system(rate_functions, params, 
                                         init_dist, unique(targets$agep))
            sum(dnorm(x = as.vector(predicted),
                      mean = targets$pct,
                      sd = targets$pct_se,
                      log = TRUE))
        })
    }
    
    # Stop here and return diagnostic information
    return(list(
        transitions = transitions,
        rate_functions = rate_functions,
        test_results = test_results,
        best_scale = best_scale,
        parameters = list(
            total = total_params,
            per_transition = n_params_per_transition,
            n_transitions = n_transitions
        )
    ))
}
calibrated <- calibrate_transitions(model, targets)






# Get posterior predictive distribution
posterior_pred <- function(calibrated, ages) {
    # Sample from posterior
    n_samples <- 100
    sample_idx <- sample(length(calibrated$weights), 
                         n_samples, 
                         prob = calibrated$weights)
    
    # Generate predictions for each sample
    predictions <- map_dfr(sample_idx, function(i) {
        pred <- simulate_system(
            calibrated$rate_functions,
            calibrated$parameters[i,],
            initial_dist = targets %>% 
                filter(agep == min(agep)) %>% 
                pull(pct),
            ages = ages
        )
        
        as_tibble(pred) %>%
            mutate(sample = i,
                   age = ages)
    })
    
    return(predictions)
}