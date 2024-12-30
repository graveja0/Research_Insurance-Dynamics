if (!exists("p_pre"))
    source(here::here("R/download-and-prepare-meps-data.r"))
if (!exists("params"))
    source(here::here("R/define-parameters.r"))
source(here::here("R/model-functions.r"))

knot_locations <- c(0,1,5, 10, 15, 18, 19, 20, 21, 26, 35, 50, 60)
knot_locations <- knot_locations[between(knot_locations,params_$age0,params_$age0+params_$horizon)]

prepare_msm_data <- function(data, knot_locations) {
    # Ensure states are numeric for msm
    states <- sort(unique(data$type))
    state_to_num <- seq_along(states)
    names(state_to_num) <- states
    
    # Create spline basis with specified knots
    age_range <- range(data$age)
    ns_basis <- ns(data$age, 
                   knots = knot_locations[-c(1, length(knot_locations))],
                   Boundary.knots = c(knot_locations[1], knot_locations[length(knot_locations)]))
    
    # Get spline basis values
    spline_matrix <- predict(ns_basis, data$age)
    colnames(spline_matrix) <- paste0("age_spline", 1:ncol(spline_matrix))
    
    # Create the model matrix for splines and scaled weight
    spline_data <- as.data.frame(spline_matrix)
    
    # Scale weight to mean 0, sd 1
    weight_scaled <- scale(data$weight)
    spline_data$weight_scaled <- as.vector(weight_scaled)
    
    # Combine with original data
    msm_data <- data %>%
        mutate(
            state = state_to_num[type]  # Convert states to numeric
        ) %>%
        bind_cols(spline_data)
    
    return(list(
        data = msm_data,
        states = states,
        state_to_num = state_to_num,
        ns_basis = ns_basis,
        knots = knot_locations,
        spline_names = c(colnames(spline_matrix), "weight_scaled"),
        weight_mean = mean(data$weight),
        weight_sd = sd(data$weight)
    ))
}

#' Predict transition probabilities for a given age and weight
# Update predict_transitions to match
predict_transitions <- function(model, age, weight = NULL) {
    Q <- get_transition_intensities(model, age, weight)
    P <- pmatrix.msm(model, t = 1, 
                     covariates = list(weight_scaled = if(!is.null(weight)) {
                         (weight - attr(model, "weight_mean")) / attr(model, "weight_sd")
                     } else {
                         0
                     }))
    
    states <- attr(model, "states")
    rownames(P) <- colnames(P) <- states
    
    return(P)
}
#' Fit multi-state model with spline-based age effects
fit_insurance_msm <- function(prepared_data, control = list(maxit = 1000, reltol = 1e-8, fnscale = 4000)) {
    # Create formula for spline terms and weight
    covariates <- prepared_data$spline_names  # Now includes weight_scaled
    covariate_formula <- paste(covariates, collapse = " + ")
    
    # Create transition matrix for initial values
    n_states <- length(prepared_data$states)
    Q <- matrix(0.1, n_states, n_states)
    diag(Q) <- -rowSums(Q)
    
    # Create a simple constraint structure
    # Format: list of numeric vectors
    n_params <- length(covariates)
    
    constraint <- list(
        # Constraint for Public (3) to Employer (1) transitions
        c(rep(1, n_params)),  # All coefficients positive
        
        # Constraint for Public (3) to OthPrivate (2) transitions
        c(rep(1, n_params)),  # All coefficients positive
        
        # Constraint for Employer (1) to Public (3) transitions
        c(rep(-1, n_params))  # All coefficients negative
    )
    
    # Fit model
    model <- msm(
        state ~ month,
        subject = patient_id,
        data = prepared_data$data,
        qmatrix = Q,
        covariates = formula(paste("~", covariate_formula)),
        control = control,
        method = "BFGS",
        constraint = constraint
    )
    
    # Add useful attributes
    attr(model, "states") <- prepared_data$states
    attr(model, "state_to_num") <- prepared_data$state_to_num
    attr(model, "knots") <- prepared_data$knots
    attr(model, "spline_names") <- prepared_data$spline_names
    
    return(model)
}

#' Plot transition probabilities with knot locations
#' @param model Fitted msm model
#' @param ages Vector of ages at which to evaluate transitions
#' @param weight Optional weight value to use for predictions
plot_transitions <- function(model, ages = seq(0, 65, by = 1)) {
    # Calculate transition probabilities at each age
    transitions <- map_dfr(ages, function(age) {
        P <- predict_transitions(model, age, weight)
        tibble(
            age = age,
            from_state = rep(rownames(P), each = ncol(P)),
            to_state = rep(colnames(P), times = nrow(P)),
            probability = as.vector(P)
        )
    })
    
    # Get knot locations
    knots <- attr(model, "knots")
    
    # Create plot
    transitions %>%
        filter(from_state != to_state) %>%  # Remove diagonal elements
        ggplot(aes(x = age, y = probability, color = to_state)) +
        geom_vline(xintercept = knots, linetype = "dashed", alpha = 0.3) +
        geom_line() +
        facet_wrap(~from_state) +
        theme_minimal() +
        labs(
            title = "Insurance Transition Probabilities by Age",
            subtitle = paste("Knots at ages:", paste(knots, collapse = ", ")),
            x = "Age",
            y = "Monthly Transition Probability",
            color = "To State"
        )
}

#' Optimize initial distribution to match targets using DEoptim
#' @param model MSM model
#' @param targets Target data frame
#' @param age_range Ages to consider in optimization
optimize_initial_dist <- function(model, targets, age_range = 0:64) {
    # Objective function for optimization
    objective <- function(x) {
        # Normalize proportions to sum to 1
        props <- x / sum(x)
        names(props) <- attr(model, "states")
        
        # Run simulation with these initial proportions
        sim_results <- simulate_lifetime(
            model = model,
            initial_dist = props,
            age_range = age_range
        )
        
        # Compare to targets
        # Calculate mean squared error at each age
        error <- sim_results %>%
            left_join(targets, by = c("age" = "agep", "type" = "type")) %>%
            mutate(sq_error = (probability - pct)^2) %>%
            summarise(mse = mean(sq_error, na.rm = TRUE)) %>%
            pull(mse)
        
        return(error)
    }
    
    # Set up DEoptim control parameters
    control <- DEoptim.control(
        NP = 100,        # Population size
        itermax = 10,   # Maximum iterations
        F = 0.8,         # Differential weight
        CR = 0.9,        # Crossover probability
        strategy = 2,    # DE strategy
        trace = TRUE     # Show progress
    )
    
    # Run optimization with DEoptim
    result <- DEoptim(
        fn = objective,
        lower = rep(0.00000000000000001, 4),  # Lower bounds
        upper = rep(0.99, 4),  # Upper bounds
        control = control
    )
    
    # Return normalized optimal proportions
    opt_props <- result$optim$bestmem / sum(result$optim$bestmem)
    names(opt_props) <- attr(model, "states")
    
    return(list(
        proportions = opt_props,
        objective_value = result$optim$bestval,
        convergence = result$convergence,
        full_result = result  # Include full optimization result for diagnostics
    ))
}



#' Get transition intensities from fitted msm model
#' @param model Fitted msm model
#' @param age Age at which to compute intensities
#' @return Matrix of transition intensities
#' Get transition intensities from fitted msm model
#' @param model Fitted msm model
#' @param age Age at which to compute intensities
#' @return Matrix of transition intensities
# Update get_transition_intensities to handle weight
get_transition_intensities <- function(model, age, weight = NULL) {
    # Get spline basis for this age
    knots <- attr(model, "knots")
    ns_basis <- ns(age, 
                   knots = knots[-c(1, length(knots))],
                   Boundary.knots = c(knots[1], knots[length(knots)]))
    
    # Create covariate vector
    spline_values <- predict(ns_basis, age)
    spline_names <- attr(model, "spline_names")
    
    # Create covariates list with spline values
    covariates <- setNames(as.list(spline_values), 
                           spline_names[!spline_names == "weight_scaled"])
    
    # Add weight if provided, otherwise use mean (0)
    covariates$weight_scaled <- if(!is.null(weight)) {
        (weight - attr(model, "weight_mean")) / attr(model, "weight_sd")
    } else {
        0
    }
    
    # Get Q matrix
    states <- attr(model, "states")
    n_states <- length(states)
    
    # Get raw qmatrix values
    q_values <- qmatrix.msm(model, covariates = covariates)
    
    # Convert to matrix form
    Q <- matrix(0, n_states, n_states)
    for(i in 1:n_states) {
        for(j in 1:n_states) {
            if(i != j) {
                Q[i,j] <- as.numeric(q_values[i,j][1])
            }
        }
    }
    
    diag(Q) <- -rowSums(Q)
    rownames(Q) <- colnames(Q) <- states
    
    return(Q)
}
#' Create annual transition probability matrix from monthly intensities
#' @param Q Monthly transition intensity matrix
#' @return Annual transition probability matrix
monthly_to_annual <- function(Q) {
    # Convert monthly to annual using matrix exponential
    P_annual <- expm(Q * 12)
    return(P_annual)
}

#' Simulate lifetime trajectories using age-varying transition matrices
#' @param model Fitted msm model
#' @param initial_dist Initial distribution across states
#' @param age_range Vector of ages to simulate
#' @return Data frame of state occupancy probabilities by age
simulate_lifetime <- function(model, 
                              initial_dist =targets %>% filter(agep==0) %>% pull(pct),
                              age_range = 1:64) {
    
    states <- attr(model, "states")
    n_states <- length(states)
    
    
    # Default initial distribution if not provided
    if(is.null(initial_dist)) {
        initial_dist <- targets %>% filter(agep==0) %>% pull(pct)
        names(initial_dist) <- states
    }
    
    # Initialize results matrix
    n_ages <- length(age_range)
    results <- matrix(0, nrow = n_ages, ncol = n_states)
    colnames(results) <- states
    
    # Set initial distribution
    current_dist <- initial_dist
    
    # Simulate forward
    for(i in seq_along(age_range)) {
        age <- age_range[i]
        
        # Store current distribution
        results[i,] <- current_dist
        
        # Get transition matrix for this age
        Q <- get_transition_intensities(model, age)
        P <- monthly_to_annual(Q)
        
        # Update distribution
        current_dist <- current_dist %*% P
    }
    
    # Convert to data frame
    results_df <- as.data.frame(results) %>%
        mutate(age = age_range) %>%
        pivot_longer(-age, 
                     names_to = "type",
                     values_to = "probability")
    
    return(results_df)
}

#' Plot simulated trajectories against calibration targets
#' @param simulated_data Output from simulate_lifetime
#' @param target_data Calibration target data
plot_comparison <- function(simulated_data, target_data) {
    # Create combined plot
    ggplot() +
        # Plot calibration targets
        geom_line(data = target_data,
                  aes(x = agep, y = pct, linetype = "Calibration Target")) +
        # Plot simulated trajectories
        geom_line(data = simulated_data,
                  aes(x = age, y = probability, linetype = "Model Prediction"),
                  color = "red") +
        facet_wrap(~type) +
        theme_minimal() +
        labs(title = "Insurance Coverage Trajectories: Model vs Calibration Targets",
             x = "Age",
             y = "Proportion",
             linetype = "Source") +
        scale_linetype_manual(values = c("solid", "dashed")) +
        theme(legend.position = "bottom")
}

#' Plot uncalibrated model transitions
#' @param model Fitted msm model
#' @param %>% Target data for comparison
plot_uncalibrated_comparison <- function(model, targets) {
    # Generate predictions using equal initial distribution
    
    observed_dist = targets %>% filter(agep==0) %>% pull(pct)
    
    # Get results from model
    uncalibrated_results <- simulate_lifetime(
        model = model,
        initial_dist = observed_dist,
        age_range = min(targets$agep):max(targets$agep)
    )
    
    # Create plot comparing uncalibrated results to targets
    ggplot() +
        # Plot calibration targets
        geom_line(data = targets,
                  aes(x = agep, y = pct, linetype = "Target")) +
        # Plot model predictions
        geom_line(data = uncalibrated_results,
                  aes(x = age, y = probability, linetype = "Uncalibrated Model"),
                  color = "blue") +
        facet_wrap(~type) +
        theme_minimal() +
        labs(
            title = "Insurance Coverage: Uncalibrated Model vs Targets",
            x = "Age",
            y = "Proportion",
            linetype = "Source"
        ) +
        scale_linetype_manual(values = c("solid", "dashed")) +
        theme(legend.position = "bottom")
}


####################
# ACS-Based Targets
####################

targets <- 
    read_rds(here::here(glue::glue("_inputs/acs-calibration-targets/acs-calibration-targets_2012.rds"))) %>% 
    mutate(type = case_when(type==1 ~ "Employer",
                            type==2 ~ "OthPrivate",
                            type==3 ~ "Public", 
                            type==4 ~ "Uninsured")) %>% 
    filter(agep >= params$age0 & agep <= params$age0 + params$horizon)

targets %>% 
    ggplot(aes(x = agep, y = pct)) + geom_line() + facet_wrap(~type) + 
    ggthemes::theme_economist_white()


targets %>% 
    ggplot(aes(x = agep, y = pct)) + 
    geom_line() + 
    geom_vline(xintercept = knot_locations, 
               linetype = "dashed", 
               color = "red", 
               alpha = 0.5) +
    facet_wrap(~type) + 
    ggthemes::theme_economist_white() +
    labs(title = "Insurance Coverage by Age",
         subtitle = "Red dashed lines show spline knot locations",
         x = "Age",
         y = "Percentage")

####################
# MEPS Input Data
####################

df =
    df_final_pre %>%
    group_by(dupersid) %>%
    mutate(month = (year - min(as.numeric(pre_lut))) * 12 + month(month)) %>%
    # Exclusion: no nonelderly medicare
    mutate(nonelderly_medicare = max(type == 5 & age < 65)) %>%
    filter(nonelderly_medicare == 0) %>%
    group_by(dupersid) %>%
    # Exclusion: 3+ months observed
    filter(n() >= 3)  %>%
    # Exclusion: age restriciton.
    filter(age >= params$age0 & age <= params$age0 + params$horizon) %>% 
    mutate(type = paste0(factor(
        type, levels = 1:4, labels = params$v_tr_names
    ))) %>%
    select(patient_id = dupersid,
           weight = longwt,
           month,
           type,
           age = age)

set.seed(123)
M = 1000
sampled_ids <- sample(unique(df$patient_id), M)
df_model <- df %>% filter(patient_id %in% sampled_ids)

# Fit model
# First fit model with more knots
model <- fit_insurance_msm(
    prepared_data = prepare_msm_data(
        data = df_model,
        knot_locations = knot_locations  # More detailed knots
    )
)
# Extract transition rates from fitted model
ages_seq <- seq(params$age0, params$age0 + params$horizon, by=0.5)
get_transition_rates <- function(model, ages) {
    # Get rates at each age
    rates_by_age <- map_dfr(ages, function(age) {
        Q <- qmatrix.msm(model, covariates = list(age = age))
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

fitted_rates <- get_transition_rates(model, ages_seq)

# Plot transition rates over age
rate_plot <- ggplot(fitted_rates, aes(x=age, y=rate, color=paste(from_state, "->", to_state))) +
    geom_line() +
    facet_wrap(~from_state, scales="free_y") +
    theme_bw() +
    labs(title="Fitted Transition Rates by Age",
         x="Age", 
         y="Rate",
         color="Transition") +
    theme(legend.position="bottom")
print(rate_plot)

# Look at model summary
print(summary(model))

# Check model convergence
print("Model convergence diagnostics:")
print(model$opt$convergence)  # 0 indicates successful convergence
print(model$opt$message)

# Print likelihood
cat("\nFinal log-likelihood:", model$minus2loglik/-2)

# Calculate AIC and BIC
n_params <- length(model$estimates)
n_obs <- nrow(df_model)
aic <- model$minus2loglik + 2*n_params
bic <- model$minus2loglik + log(n_obs)*n_params

cat("\nModel fit statistics:")
cat("\nAIC:", aic)
cat("\nBIC:", bic)

# Plot observed vs expected prevalence
prevalence_plot <- plot.msm(model, legend.pos=c(8,1))
print(prevalence_plot)


# Plot uncalibrated results
uncalibrated_plot <- plot_uncalibrated_comparison(model, targets)
print(uncalibrated_plot)

