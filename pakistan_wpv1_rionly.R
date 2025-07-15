# SEIR model for WPV1 in Pakistan with:
#   province-specific compartments, mixing
#   seasonal forcing
#   routine immunization
#   daily births

#### Libraries ####
library(deSolve)
library(tidyverse)
library(adaptivetau)
library(matrixStats)

#### Model Setup ####
# Define provinces
provinces <- c("Punjab", "Sindh", "KPK", "KPK_Sub", "Balochistan", "Balochistan_Sub")

# Seasonal forcing function
seasonal_beta <- function(beta, time, amp, peak_day) {
  beta * (1 + amp * cos(2 * pi * (time - peak_day) / 365))
}

# Parameters
params <- list(
  beta = c(Punjab = 0.57, Sindh = 0.57, KPK = 0.57, KPK_Sub = 0.57, Balochistan = 0.57, Balochistan_Sub = 0.57), # Basic reproductive numbers fixed at ~8
  sigma = 1/7,  # Latent period ~7 days
  gamma = 1/14, # Infectious period ~14 days
  vac_frac = c(Punjab = 0.83, Sindh = 0.71, KPK = 0.67, KPK_Sub = 0.3, Balochistan = 0.52, Balochistan_Sub = 0.3), # fraction protected by routine immunization (assumed coverage * take). IDM immunity mapper dpt1 * per-dose take rate 
  seasonal_amp = 0.15, # seasonal amplitude
  peak_day = 180, # June 30
  init_pop = c(Punjab = 17000000, Sindh = 8000000, KPK = 3800000, KPK_Sub = 250000, Balochistan = 1400000, Balochistan_Sub = 400000), # Initial population for each province. Allocate the estimated 5% of pakistan pop into sub-pops proportionally into KPK and Balochistan
  imm_start = c(Punjab = 0.9734, Sindh = 0.9612, KPK = 0.9544, KPK_Sub = 0.9, Balochistan = 0.9245, Balochistan_Sub = 0.9), # Immunization coverage at Jan 2026 per IDM Immunity Mapper. Set sub populations to have R=1 at start
  afp_frac = 1/200
)

# Add daily number of births by each province. Also assumed to be the daily exit rate as populations age-out.
national_daily_births <- 19176
params$daily_births <- c(params$init_pop / sum(params$init_pop) * national_daily_births)

# Add a mixing matrix representing differential mixing within and between provinces
params$mixing_matrix <- matrix(
  # Pun   Sin   KPK  KPK_Sub Bal   Bal_Sub
  c(0.95, 0.01, 0.01, 0.01,  0.01, 0.01,   # Pun
    0.01, 0.95, 0.01, 0.01,  0.01, 0.01,   # Sin
    0.01, 0.01, 0.91, 0.05,  0.01, 0.01,   # KPK
    0.01, 0.01, 0.05, 0.91,  0.01, 0.01,   # KPK_Sub
    0.01, 0.01, 0.01, 0.01,  0.91, 0.05,   # Bal
    0.01, 0.01, 0.01, 0.01,  0.05, 0.91),  # Bal_Sub
  nrow = 6, byrow = TRUE,
  dimnames = list(provinces, provinces)
)

# Initial state variables (Susceptible, Exposed, Infectious, Recovered)
# Population is for children under 15 years of age
init_state <- c(
  S_Punjab = as.numeric(params$init_pop["Punjab"] * (1-params$imm_start["Punjab"])), E_Punjab = 0, I_Punjab = 140, R_Punjab = as.numeric(params$init_pop["Punjab"] * params$imm_start["Punjab"]),AFP_Punjab = 0,
  S_Sindh = as.numeric(params$init_pop["Sindh"] * (1-params$imm_start["Sindh"])), E_Sindh = 0, I_Sindh = 140, R_Sindh = as.numeric(params$init_pop["Sindh"] * params$imm_start["Sindh"]), AFP_Sindh = 0,
  S_KPK = as.numeric(params$init_pop["KPK"] * (1-params$imm_start["KPK"])), E_KPK = 0, I_KPK = 70, R_KPK = as.numeric(params$init_pop["KPK"] * params$imm_start["KPK"]), AFP_KPK = 0,
  S_KPK_Sub = as.numeric(params$init_pop["KPK_Sub"] * (1-params$imm_start["KPK_Sub"])), E_KPK_Sub = 0, I_KPK_Sub = 70, R_KPK_Sub = as.numeric(params$init_pop["KPK_Sub"] * params$imm_start["KPK_Sub"]), AFP_KPK_Sub = 0,
  S_Balochistan = as.numeric(params$init_pop["Balochistan"] * (1-params$imm_start["Balochistan"])), E_Balochistan = 0, I_Balochistan = 70, R_Balochistan = as.numeric(params$init_pop["Balochistan"] * params$imm_start["Balochistan"]), AFP_Balochistan = 0,
  S_Balochistan_Sub = as.numeric(params$init_pop["Balochistan_Sub"] * (1-params$imm_start["Balochistan_Sub"])), E_Balochistan_Sub = 0, I_Balochistan_Sub = 70, R_Balochistan_Sub = as.numeric(params$init_pop["Balochistan_Sub"] * params$imm_start["Balochistan_Sub"]), AFP_Balochistan_Sub = 0
)

#### SEIR Model Function ####
seir_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    derivatives <- c()
    
    for(i in 1:length(provinces)) {
      prov <- provinces[i]
      S <- state[paste0("S_", prov)]
      E <- state[paste0("E_", prov)]
      I <- state[paste0("I_", prov)]
      R <- state[paste0("R_", prov)]
      N <- S + E + I + R
      AFP <- state[paste("AFP_", prov)]
      
      beta_t <- seasonal_beta(beta[prov], time, seasonal_amp, peak_day)
      
      # Calculate force of infection with differential mixing
      lambda <- 0
      for(j in 1:length(provinces)){
        other_prov <- provinces[j]
        I_other <- state[paste0("I_", other_prov)]
        N_other <- sum(state[paste0(c("S_", "E_", "I_", "R_"), other_prov)])
        lambda <- lambda + mixing_matrix[prov, other_prov] * beta_t * (I_other / N_other)
      }
      
      dS <- daily_births[prov]*(1-vac_frac[prov]) - lambda * S - daily_births[prov] * S/N
      dE <- lambda * S - sigma * E - daily_births[prov] * E/N
      dI <- sigma * E - gamma * I - daily_births[prov] * I/N
      dR <- daily_births[prov]*vac_frac[prov] + gamma * I - daily_births[prov] * R/N
      dAFP <- params$afp_frac * lambda * S
      
      derivatives <- c(derivatives, dS, dE, dI, dR, dAFP)
    }
    
    names(derivatives) <- names(state)
    return(list(derivatives))
  })
}

#### Simulation ####
sim_time <- seq(0, 365 * 10, by = 1)

# Run the simulation
output <- ode(y = init_state, 
              times = sim_time, 
              func = seir_model, 
              parms = params)

output_df <- as.data.frame(output)

# Add N for each province
output_df <- output_df %>%
  mutate(N_Punjab = S_Punjab + E_Punjab + I_Punjab + R_Punjab,
         N_Sindh = S_Sindh + E_Sindh + I_Sindh + R_Sindh,
         N_KPK = S_KPK + E_KPK + I_KPK + R_KPK,
         N_KPK_Sub = S_KPK_Sub + E_KPK_Sub + I_KPK_Sub + R_KPK_Sub,
         N_Balochistan = S_Balochistan + E_Balochistan + I_Balochistan + R_Balochistan,
         N_Balochistan_Sub = S_Balochistan_Sub + E_Balochistan_Sub + I_Balochistan_Sub + R_Balochistan_Sub)

# Calculate incident AFP cases
output_df <- output_df %>%
  mutate(AFP_incident_Punjab = c(0, diff(AFP_Punjab)),
         AFP_incident_Sindh = c(0, diff(AFP_Sindh)),
         AFP_incident_KPK = c(0, diff(AFP_KPK)),
         AFP_incident_KPK_Sub = c(0, diff(AFP_KPK_Sub)),
         AFP_incident_Balochistan = c(0, diff(AFP_Balochistan)),
         AFP_incident_Balochistan_Sub = c(0, diff(AFP_Balochistan_Sub)))

#### Plot results ####
par(mfrow=c(3,2))

for(prov in provinces){
  plot(output_df$time / 365 + 2026, output_df[[paste0("S_", prov)]], type = "l", col = "blue",
       xlab = "Year", ylab = "Susceptible Pop", main = paste("Susceptibles -", prov))
}

for(prov in provinces){
  plot(output_df$time / 365 + 2026, output_df[[paste0("S_", prov)]] / output_df[[paste0("N_", prov)]], type = "l", col = "green",
       xlab = "Year", ylab = "Proportion", main = paste("Proportion Susceptible -", prov), ylim = c(0, 0.2))
  lines(x = output_df$time / 365 + 2026, y = rep(1/(params$beta[[paste0(prov)]]/params$gamma), times = length(sim_time)), col = "black", lty = "dashed")
}

for(prov in provinces){
  plot(output_df$time / 365 + 2026, output_df[[paste0("I_", prov)]], type = "l", col = "red",
       xlab = "Year", ylab = "Infections", main = paste("WPV1 Infections -", prov))
}

for(prov in provinces){
  plot(output_df$time / 365 + 2026, output_df[[paste0("AFP_", prov)]], type = "l", col = "red",
       xlab = "Year", ylab = "Cumulative AFP Cases", main = paste("AFP Cases -", prov))
}

for(prov in provinces){
  plot(output_df$time / 365 + 2026, output_df[[paste0("AFP_", prov)]] / output_df[[paste0("N_", prov)]] * 100, type = "l", col = "red",
       xlab = "Year", ylab = "Attack Rate (%)", main = paste("AFP Attack Rate -", prov))
}

for(prov in provinces){
  plot(output_df$time / 365 + 2026, output_df[[paste0("AFP_incident_", prov)]], type = "l", col = "red",
       xlab = "Year", ylab = "Daily AFP Cases", main = paste("AFP Cases -", prov))
}

# plot monthly AFP cases 
output_df$month <- floor(output_df$time / 30) + 1
output_df_monthly <- output_df %>%
  group_by(month) %>%
  summarise(across(starts_with("AFP_incident_"), sum, na.rm = TRUE))
for(prov in provinces){
  plot(output_df_monthly$month/12 + 2026, output_df_monthly[[paste0("AFP_incident_", prov)]], type = "l", col = "red",
       xlab = "Year", ylab = "Monthly AFP Cases", main = paste("AFP Cases -", prov))
}

# Plot total AFP cases across all provinces
par(mfrow=c(1,1))
total_afp_cases <- rowSums(output_df_monthly[, paste0("AFP_incident_", provinces)])
plot(output_df_monthly$month/12 + 2026, total_afp_cases, type = "l", col = "red",
     xlab = "Year", ylab = "Total Monthly AFP Cases", main = "Total AFP Cases Across Pakistan")

# Estimate endemic equilibrium annual case count
mean(tail(total_afp_cases,24)) * 12 # Monthly to annual

#### Explore inputs ####
# Basic reproductive number
params$beta / params$gamma # by province
sum((params$beta / params$gamma) * params$init_pop/sum(params$init_pop)) # weighted sum

# Effective reproductive number at start of simulation
params$beta / params$gamma * (1-params$imm_start)

# Effective reproductive number at year 1 (absent infections)
S_year1 = init_state[paste0("S_", provinces)] + params$daily_births * 365 * (1-params$vac_frac)
R_year1 = init_state[paste0("R_", provinces)] + params$daily_births * 365 * params$vac_frac
S_year1 / (S_year1 + R_year1) * params$beta / params$gamma

# Effective reproductive number at year 2 (absent infections)
S_year2 = init_state[paste0("S_", provinces)] + params$daily_births * 365*2 * (1-params$vac_frac)
R_year2 = init_state[paste0("R_", provinces)] + params$daily_births * 365*2 * params$vac_frac
S_year2 / (S_year2 + R_year2) * params$beta / params$gamma

# Effective reproductive number at year 5 (absent infections)
S_year5 = init_state[paste0("S_", provinces)] + params$daily_births * 365 * 5 * (1-params$vac_frac)
R_year5 = init_state[paste0("R_", provinces)] + params$daily_births * 365 * 5 * params$vac_frac
S_year5 / (S_year5 + R_year5) * params$beta / params$gamma

# Initialize infections, given 74 AFP cases in 2024, case:infection ratio, and applying evenly across provinces
74 * 1/params$afp_frac / (365 * params$gamma) #estimated number of infections on an average 2024 day in Pakistan

# RI estimates
  # Kid Risk models per-bOPV-dose take rate at 0.48 and national RI coverage at 0.75 (with undervaccinated subpopulations at 0.4). 
  # Four bOPV doses in Pakistan schedule.
(1-(1-0.48)^4) # protection for those on-schedule
(1-(1-0.48)^4) * .75 # average national protection
(1-(1-0.48)^4) * .4 # protection for undervaccinated subpopulations
  # ICL endemic eq analysis used 84% OPV immune.
sum(params$vac_frac * params$init_pop) / sum(params$init_pop) # average national protection by province

# Plot seasonal beta function
seasonal_times <- seq(0, 365, by = 1)
seasonal_betas <- sapply(seasonal_times, function(t) {
  seasonal_beta(params$beta["KPK"] / params$gamma, t, params$seasonal_amp, params$peak_day)
})
plot(seasonal_times, seasonal_betas, type = "l", col = "blue",
     xlab = "Day of Year", ylab = "Seasonal R_0",
     main = "Seasonal Forcing of Beta in KPK")

# Plot mixing matrix
heatmap(params$mixing_matrix, Rowv = NA, Colv = NA, 
        scale = "none", col = heat.colors(256), 
        main = "Mixing Matrix for Provinces",
        xlab = "Provinces", ylab = "Provinces")

#### Stochastic simulation ####
# Convert SEIR model into stochastic simulation using Gillespie algorithm
library(adaptivetau)

# State vector
state_names <- c()
for(prov in provinces){
  state_names <- c(state_names, paste0(c("S_", "E_", "I_", "R_", "AFP_"), prov))
}

initial_state <- setNames(as.integer(init_state), state_names)

# Define transitions
transitions <- list()
for(i in provinces){
  # Infection
  transitions[[paste0("infection_", i)]] <- setNames(c(-1,1,0,0,0), paste0(c("S_", "E_", "I_", "R_", "AFP_"), i))
  # Progression E to I
  transitions[[paste0("progression_", i)]] <- setNames(c(0,-1,1,0,0), paste0(c("S_", "E_", "I_", "R_", "AFP_"), i))
  # Recovery
  transitions[[paste0("recovery_", i)]] <- setNames(c(0,0,-1,1,0), paste0(c("S_", "E_", "I_", "R_", "AFP_"), i))
  # AFP case counting
  transitions[[paste0("AFP_case_", i)]] <- setNames(c(0,0,0,0,1), paste0(c("S_", "E_", "I_", "R_", "AFP_"), i))
  # Birth and Vaccination
  transitions[[paste0("birth_vaccinated_", i)]] <- setNames(c(0,0,0,1,0), paste0(c("S_", "E_", "I_", "R_", "AFP_"), i))
  transitions[[paste0("birth_unvaccinated_", i)]] <- setNames(c(1,0,0,0,0), paste0(c("S_", "E_", "I_", "R_", "AFP_"), i))
  # Aging out (equal rates for simplicity)
  for(comp in c("S_", "E_", "I_", "R_")){
    transitions[[paste0("aging_", comp, i)]] <- setNames(c(0,0,0,0,0), paste0(c("S_", "E_", "I_", "R_", "AFP_"), i))
    transitions[[paste0("aging_", comp, i)]][[paste0(comp,i)]] <- -1
  }
}

# Rate function
rate_func <- function(x, params, t){
  rates <- numeric()
  for(i in provinces){
    N <- sum(x[paste0(c("S_","E_","I_","R_"), i)])
    beta_t <- seasonal_beta(params$beta[i], t, params$seasonal_amp, params$peak_day)
    lambda <- 0
    for(j in provinces){
      I_other <- x[paste0("I_", j)]
      N_other <- sum(x[paste0(c("S_","E_","I_","R_"), j)])
      lambda <- lambda + params$mixing_matrix[i,j] * beta_t * (I_other/N_other)
    }
    
    rates <- c(rates,
               lambda * x[paste0("S_", i)],  # infection
               params$sigma * x[paste0("E_", i)], # progression
               params$gamma * x[paste0("I_", i)], # recovery
               params$afp_frac * lambda * x[paste0("S_", i)], # AFP cases
               params$daily_births[i]*params$vac_frac[i], # vaccinated birth
               params$daily_births[i]*(1-params$vac_frac[i]), # unvaccinated birth
               params$daily_births[i]*x[paste0("S_",i)]/N, # aging S
               params$daily_births[i]*x[paste0("E_",i)]/N, # aging E
               params$daily_births[i]*x[paste0("I_",i)]/N, # aging I
               params$daily_births[i]*x[paste0("R_",i)]/N # aging R
    )
  }
  return(rates)
}

# Simulate
sim_result <- ssa.adaptivetau(init.values = initial_state, 
                              transitions = transitions, 
                              rateFunc = rate_func,
                              params = params,
                              tf = 365*10)
sim_result <- as.data.frame(sim_result)

# Create field for incident AFP cases
for(prov in provinces){
  sim_result[,paste0("AFP_incident_", prov)] <- c(0, diff(sim_result[,paste0("AFP_", prov)]))
}

# Plot results for infections
par(mfrow=c(3,2))
for(prov in provinces){
  plot(sim_result[,"time"]/365+2026, sim_result[,paste0("I_",prov)], type="l", col="red",
       xlab="Year", ylab="Infections", main=paste("Stochastic WPV1 Infections -",prov))
}

# Repeat simulations
library(matrixStats)

n_runs <- 50
time_points <- seq(0, 365*10, by = 1)
simulation_results_I <- array(NA, dim = c(length(time_points), length(provinces), n_runs), 
                            dimnames = list(NULL, provinces, NULL))
simulation_results_AFP <- array(NA, dim = c(length(time_points), length(provinces), n_runs), 
                            dimnames = list(NULL, provinces, NULL))

for(run in 1:n_runs){
  sim_result <- ssa.adaptivetau(init.values = initial_state, 
                                transitions = transitions, 
                                rateFunc = rate_func,
                                params = params,
                                tf = 365*10)
  
  for(i in 1:length(provinces)){
    prov <- provinces[i]
    inf_res <- approx(sim_result[,"time"], sim_result[,paste0("I_", prov)], xout = time_points)$y
    simulation_results_I[,prov,run] <- inf_res
    
    afp_res <- approx(sim_result[,"time"], sim_result[,paste0("AFP_", prov)], xout = time_points)$y
    simulation_results_AFP[,prov,run] <- afp_res
  }
  cat(".")
}

simulation_matrix_I <- simplify2array(simulation_results_I)
simulation_matrix_AFP <- simplify2array(simulation_results_AFP)

# Convert cumulative AFP to monthly incident AFP events
simulation_matrix_AFP_incident <- array(NA, dim = dim(simulation_matrix_AFP))
simulation_matrix_AFP_incident_monthly <- array(NA, dim = c(1 + length(time_points) %/% 30, length(provinces), n_runs),
                                             dimnames = list(NULL, provinces, NULL))
for(i in 1:length(provinces)){
  for(run in 1:n_runs){
    simulation_matrix_AFP_incident[,i,run] <- c(0, diff(simulation_matrix_AFP[,i,run]))
    # Convert daily into monthly
    simulation_matrix_AFP_incident_monthly[,i,run] <- tapply(simulation_matrix_AFP_incident[,i,run], 
                                                     (1:length(time_points) - 1) %/% 30, sum)
  }
}

# Plot I, by province
par(mfrow=c(3,2))
for(i in 1:length(provinces)){
  prov <- provinces[i]
  mat <- simulation_matrix_I[, i, ]
  plot(time_points/365+2026, mat[,1], type="l", col=rgb(1,0,0,0.2),
       xlab="Year", ylab="Infections", main=paste("Stochastic WPV1 Infections -", prov), ylim=range(mat))
  for(run in 2:n_runs){
    lines(time_points/365+2026, mat[,run], col=rgb(1,0,0,0.2))
  }
  median_trace <- rowMedians(mat)
  lines(time_points/365+2026, median_trace, col="black", lwd=2)
}

# Plot AFP, by province, by month
par(mfrow=c(3,2))
for(i in 1:length(provinces)){
  prov <- provinces[i]
  mat <- simulation_matrix_AFP_incident_monthly[, i, ]
  plot(1:nrow(mat) / 12 + 2026, mat[,1], type="l", col=rgb(1,0,0,0.2),
       xlab="Year", ylab="Monthly AFP Cases", main=paste("Stochastic WPV1 AFP -", prov), ylim=range(mat))
  for(run in 2:n_runs){
    lines(1:nrow(mat) / 12 + 2026, mat[,run], col=rgb(1,0,0,0.2))
  }
  median_trace <- rowMedians(mat)
  lines(1:nrow(mat) / 12 + 2026, median_trace, col="black", lwd=2)
}

# Plot monthly AFP across country
simulation_matrix_AFP_incident_monthly_country <- array(NA, dim = c(1 + length(time_points) %/% 30, n_runs),
                                                dimnames = list(NULL, NULL))
for(run in 1:n_runs){
  simulation_matrix_AFP_incident_monthly_country[,run] <- rowSums(simulation_matrix_AFP_incident_monthly[,,run])
}

# Plot monthly AFP across country
par(mfrow=c(1,1))
plot(1:nrow(simulation_matrix_AFP_incident_monthly_country) / 12 + 2026, 
     simulation_matrix_AFP_incident_monthly_country[,1], type="l", col=rgb(1,0,0,0.2),
     xlab="Year", ylab="Monthly AFP Cases", main="Stochastic WPV1 AFP - Pakistan", ylim=range(simulation_matrix_AFP_incident_monthly_country))
for(run in 2:n_runs){
  lines(1:nrow(simulation_matrix_AFP_incident_monthly_country) / 12 + 2026, 
        simulation_matrix_AFP_incident_monthly_country[,run], col=rgb(1,0,0,0.2))
}
median_trace_country <- rowMedians(simulation_matrix_AFP_incident_monthly_country)
lines(1:nrow(simulation_matrix_AFP_incident_monthly_country) / 12 + 2026, median_trace_country, col="black", lwd=2)

# Estimate endemic equilibrium annual case count
mean(tail(median_trace_country,36)) * 12 # Monthly to annual

# Proportion of simulations reaching die-out within a year
mean(apply(simulation_matrix_I[365, , ], 2, function(x) all(x < 1)))

