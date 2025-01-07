## Function to fit multiple functional forms to survival data
# without GenGamma: Reaction DARTH workgroup on GenGamma: "We removed it because the implementation in R is a bit buggy. 
# Generalized gamma density function is defined using an incomplete gamma in the denominator and when the values of the parameters 
# approach the extreme values (e.g. close to 0) then the function fails"

#Function to determine which health state occurs
samplev <- function(m.Probs, m) {
  # Arguments
  # m.Probs: matrix with probabilities (n.i * n.s)
  # m:       number of states than need to be sampled per individual  
  # Return
  # ran:    n.i x m matrix filled with sampled health state(s) per individual
  
  d <- dim(m.Probs)  # dimensions of the matrix filled with the multinomical probabilities for the health states 
  n <- d[1]          # first dimension - number of rows (number of individuals to sample for)
  k <- d[2]          # second dimension - number of columns (number of health states considered)
  lev <- dimnames(m.Probs)[[2]]  # extract the names of the health states considered for sampling
  if (!length(lev))  # in case names for the health states are missing, use numbers to specify the health states
    lev <- 1:k       # create a sequence from 1:k (number of health states considered)
  # create a matrix 
  ran <- matrix(lev[1], ncol = m, nrow = n) # create the matrix ran, filled with the first health state of the levels 
  U <- t(m.Probs)    # transposed m.Probs matrix n.i x n.s --> n.s x n.i 
  
  for(i in 2:k) {    # start loop, from the 2nd health states
    U[i, ] <- U[i, ] + U[i - 1, ] # start summing the probabilities of the different health states per individual 
  }
  if (any((U[k, ] - 1) > 1e-05))  # sum of all probs per individual - 1 should be 0 (use 1e-05 for rounding issues), else print the error statement
    stop("error in multinom: probabilities do not sum to 1")
  
  for (j in 1:m) {   # start loop of the state that needs to be sampled (m)
    un <- rep(runif(n), rep(k, n))       # sample from a uniform distribution of length n*k
    ran[, j] <- lev[1 + colSums(un > U)] # store the health state at the jth column of the U matrix
  }
  ran # return the new health state per individual n.i x m
} # close the function 

#------------------------------------------------------------------------------#
#### R function to Function to plot health state trace                      ####
####                                                                        ####
#------------------------------------------------------------------------------#

#
plot_m_TR <- function(m_M) {
  # plot the distribution of the population across health states over time (trace)
  # count the number of individuals in each health state at each cycle
  m_TR <- t(apply(m_M, 2, function(x) table(factor(x, levels = v_n, ordered = TRUE)))) 
  m_TR <- m_TR / n_i                                       # calculate the proportion of individuals 
  colnames(m_TR) <- v_n                                    # name the rows of the matrix
  rownames(m_TR) <- paste("Cycle", 0:n_t, sep = " ")       # name the columns of the matrix
  # Plot trace of first health state
  matplot(m_TR, type = "l", main = "Health state trace", col= 1:n_s,
          ylim = c(0, 1), ylab = "Proportion of cohort", xlab = "Cycle")
  legend("topright", v_n, col = 1:n_s,    # add a legend to current plot
         lty = rep(1, 3), bty = "n", cex = 0.65)
  
}


#------------------------------------------------------------------------------#
#### R function to extract the parameters of a beta distribution            ####
####                from mean and st. deviation                             ####
#------------------------------------------------------------------------------#
#' @param m mean 
#' @param s standard deviation
#' 
betaPar <- function(m, s) 
{
  a <- m * ((m * (1 - m) / s ^ 2) - 1)
  b <- (1 - m) * ((m * (1 - m) / s ^ 2) - 1)
  list(a = a, b = b)
}



#------------------------------------------------------------------------------#
#### R function to extract the parameters of a gamma distribution           ####
####                   from mean and st. deviation                          ####
#------------------------------------------------------------------------------#
#' @param m mean 
#' @param s standard deviation
#' 
gammaPar <- function(m, s) {   
  # m: mean  
  # s: standard deviation 
  shape <- m ^ 2 / s ^ 2
  scale <- s ^ 2 / m
  list(shape = shape, scale = scale)
}


#------------------------------------------------------------------------------#
#### R function to treatment lines over time                                ####
####                                                                        ####
#------------------------------------------------------------------------------#

#Plot 
plot_trace_line <- function(m_L) 
{
  m_TR <- t(apply(m_L, 2, function(x) table(factor(x, levels = c(1,2,3,4,5, 444, 555, 999), 
                                                   ordered = TRUE))))
  m_TR <- m_TR/n_i
  colnames(m_TR) <- v_names_lines
  rownames(m_TR) <- paste("Cycle", 0:n_cycles, sep = " ")
  plot(0:n_cycles, m_TR[, 1], type = "l", main = "Line membership trace", 
       ylim = c(0, 1), ylab = "Proportion of cohort", xlab = "Cycle")
  for (n_states in 2:length(v_names_lines)) {
    lines(0:n_cycles, m_TR[, n_states], col = n_states)
  }
  legend("topright", v_names_lines, col = 1:length(v_names_lines), 
         lty = rep(1, length(v_names_lines)), bty = "n", cex = 0.65)
}

plot_markov <- function (m_M) 
{
  m_TR <- t(apply(m_M, 2, function(x) table(factor(x, levels = v_names_states, 
                                                   ordered = TRUE))))
  m_TR <- m_TR/n_i
  colnames(m_TR) <- v_names_states
  rownames(m_TR) <- paste("Cycle", 0:(ncol(m_M) - 1), sep = " ")
  plot(0:(ncol(m_M) - 1), m_TR[, 1], type = "l", main = "Health state trace", 
       ylim = c(0, 1), ylab = "Proportion of cohort", xlab = "Cycle")
  for (n_states in 2:length(v_names_states)) {
    lines(0:(ncol(m_M) - 1), m_TR[, n_states], col = n_states)
  }
  legend("topright", v_names_states, col = 1:length(v_names_states), 
         lty = rep(1, length(v_names_states)), bty = "n", cex = 0.65)
}


#------------------------------------------------------------------------------#
#### R function to calculate life expectancy                                ####
#### returns life table from mortality rate input                           ####
#------------------------------------------------------------------------------#

f_le <- function(r_d){       
  
  #r_d <- mortality rate from 0 to 100
  #Age <- the age for which the life expectancy has to be returned
  
  #definitions
  # mx = mortality rate at age x
  # qx = probability to die within a year at age x
  # lx = persons alive at beginning of interval
  # dx = persons that have died in the interval
  # Lx = life years lived between successive ages
  # Tx = cumulative years lived after age x
  # ex = life expectancy at age x
  
  df_r_mort <- data.frame(mx = r_d) 
  df_r_mort$Age <- seq(from = 0, to = (length(r_d)-1)) #add age up to a hundred
  df_r_mort$qx <- 1-exp(-df_r_mort$mx) #change rate to prob
  df_r_mort$lx <- 1000  #temporary, will be changed with loop hereafter
  
  #loop to get lx, i.e. reduce lx with the number of deaths from the previous period
  for(i in 1:(nrow(df_r_mort)-1)){
    df_r_mort$lx[i+1] <- df_r_mort$lx[i]-(df_r_mort$qx[i]*df_r_mort$lx[i])  
  }
  
  df_r_mort$dx <- df_r_mort$qx*df_r_mort$lx #people that died
  
  
  df_r_mort$Lx_temp_1 <- df_r_mort$lx[2] + df_r_mort$dx[1]*0.2  #adjust first year with different half cycle correction due to infant mortality
  
  for(i in 2:(nrow(df_r_mort))-1){
    df_r_mort$Lx_temp[i] <- df_r_mort$lx[i+1] + df_r_mort$dx[i]*0.5 
    
  }
  
  
  df_r_mort$Lx_temp[1] <- df_r_mort$lx[2] + df_r_mort$dx[1]*0.2  #replace first half cycle correction value
  df_r_mort$Lx_temp[nrow(df_r_mort)] <- df_r_mort$dx[nrow(df_r_mort)]*0.5 #replace last value
  df_r_mort$Lx <- df_r_mort$Lx_temp
  df_r_mort <- df_r_mort[order(-df_r_mort$Age),] #reorder for the loop below
  
  df_r_mort$Tx[1] <- df_r_mort$Lx[1] 
  
  for(i in 1:(nrow(df_r_mort)-1)){
    df_r_mort$Tx[i+1] <- df_r_mort$Tx[i]+df_r_mort$Lx[i+1]
    
  }
  
  df_r_mort <- df_r_mort[order(df_r_mort$Age),] #reorder for the loop below
  
  df_r_mort$ex <- df_r_mort$Tx/df_r_mort$lx
  
  return(df_r_mort)
}



#------------------------------------------------------------------------------#
#### R function to extracts survival probabilities for each combination     ####
#### of covariates from the norm.boot sampled survival functions            ####
#------------------------------------------------------------------------------#

#PSA function
#MV: this function helps extract parameters from the survival objects for the PSA 
# f_extract_tp_PSA <- function(newdata, survival_object, simulated_object, n_sim, distribution){ 
#   for(i in 1:nrow(newdata)){
#     
#     df_tp_FL1_PSA <- data.frame(time = times)
#     df_tp_FL1_PSA[,2:(1+ncol(newdata))] <- newdata[i,]
#     
#     
#     for(j in 1:n_sim){
#       #For each individual i, with j (= n_sim) parameters for the PSA, replace these parameters in the survival object and extract transition probabilities
#       
#       if(distribution == "gamma"){
#         survival_object$res.t[1] <- log(simulated_object[[i]][j,1]) #shape, note distribution specific
#         survival_object$res.t[2] <- log(simulated_object[[i]][j,2]) #scale, note distribution specific
#       }
#       
#       if(distribution == "lnorm"){
#         survival_object$res.t[1] <- (simulated_object[[i]][j,1]) #meanlog, note distribution specific
#         survival_object$res.t[2] <- log(simulated_object[[i]][j,2]) #sdlog, note distribution specific
#       }
#       
#       if(distribution == "gompertz"){
#         survival_object$res.t[1] <- (simulated_object[[i]][j,1]) #shape, note distribution specific
#         survival_object$res.t[2] <- log(simulated_object[[i]][j,2]) #rate, note distribution specific
#       }
#       
#       #if(distribution == "exponential"){ #MV: does not work yet,  not needed either
#       #  survival_object$res.t[1] <- log(simulated_object[[i]][j,1]) #rate, note distribution specific
#       #}
#       
#       cumhaz_covs1 <- summary(survival_object, t = times, newdata = newdata[i,], type = "survival", tidy = TRUE)
#       df_tp_FL1_PSA[,(j+ncol(newdata)+1)] <- c(NA,trans_prob(cumhaz_covs1$est))
#       
#     }
#     #store for individual i the j PSA parameters in list number i
#     l_PSA_input_FL1[[i]] <- df_tp_FL1_PSA
#   }
#   df_tp_FL1_PSA <- as.data.frame(rbindlist(l_PSA_input_FL1)) #this gives the long format result of the list (i.e. times*length(all_covariates) long, with a column for each run of the PSA)
#   df_tp_FL1_PSA <- df_tp_FL1_PSA[,-c(1:7)]
#   return(df_tp_FL1_PSA)
# }

# f_extract_tp_PSA <- function(newdata, survival_object, simulated_object, n_sim, distribution, times) {
#   # Pre-allocate the result list
#   l_PSA_input_FL1 <- vector("list", nrow(newdata))
#   
#   # Create a lookup for distribution-specific operations
#   dist_ops <- list(
#     gamma = function(x) c(log(x[1]), log(x[2])),
#     lnorm = function(x) c(x[1], log(x[2])),
#     gompertz = function(x) c(x[1], log(x[2]))
#   )
#   
#   # Use lapply instead of for loop
#   l_PSA_input_FL1 <- lapply(seq_len(nrow(newdata)), function(i) {
#     df_tp_FL1_PSA <- data.frame(time = times)
#     df_tp_FL1_PSA[, 2:(1 + ncol(newdata))] <- newdata[i, ]
#     
#     # Pre-allocate the results matrix
#     results_matrix <- matrix(NA, nrow = length(times), ncol = n_sim)
#     
#     # Use vapply instead of for loop
#     results_matrix <- vapply(seq_len(n_sim), function(j) {
#       # Apply distribution-specific operations
#       survival_object$res.t[1:2] <- dist_ops[[distribution]](simulated_object[[i]][j, 1:2])
#       
#       cumhaz_covs1 <- summary(survival_object, t = times, newdata = newdata[i, ], type = "survival", tidy = TRUE)
#       c(NA, trans_prob(cumhaz_covs1$est))
#     }, numeric(length(times)))
#     
#     # Combine results
#     cbind(df_tp_FL1_PSA, results_matrix)
#   })
#   
#   # Combine results and remove unnecessary columns
#   df_tp_FL1_PSA <- do.call(rbind, l_PSA_input_FL1)
#   df_tp_FL1_PSA[, -(1:7)]
# }

#MV: faster parallelized function improves speed 10x
f_extract_tp_PSA <- function(newdata, survival_object, simulated_object, n_sim, distribution, times, n_cores = detectCores() - 1) {
  # Create a lookup for distribution-specific operations
  dist_ops <- list(
    gamma = function(x) c(log(x[1]), log(x[2])),
    lnorm = function(x) c(x[1], log(x[2])),
    gompertz = function(x) c(x[1], log(x[2]))
  )
  
  # Create cluster
  cl <- makeCluster(n_cores)
  
  # Define the worker function
  worker_function <- function(i, newdata, survival_object, simulated_object, n_sim, distribution, times, dist_ops) {
    df_tp_FL1_PSA <- data.frame(time = times)
    df_tp_FL1_PSA[, 2:(1 + ncol(newdata))] <- newdata[i, ]
    
    # Pre-allocate the results matrix
    results_matrix <- matrix(NA, nrow = length(times), ncol = n_sim)
    
    for (j in seq_len(n_sim)) {
      # Apply distribution-specific operations
      survival_object$res.t[1:2] <- dist_ops[[distribution]](simulated_object[[i]][j, 1:2])
      
      cumhaz_covs1 <- summary(survival_object, t = times, newdata = newdata[i, ], type = "survival", tidy = TRUE)
      results_matrix[, j] <- c(NA, trans_prob(cumhaz_covs1$est))
    }
    
    # Combine results
    cbind(df_tp_FL1_PSA, results_matrix)
  }
  
  # Export necessary functions to all workers
  clusterExport(cl, c("trans_prob"))
  
  # Parallel computation
  l_PSA_input_FL1 <- parLapply(cl, seq_len(nrow(newdata)), worker_function, 
                               newdata = newdata, 
                               survival_object = survival_object, 
                               simulated_object = simulated_object, 
                               n_sim = n_sim, 
                               distribution = distribution, 
                               times = times, 
                               dist_ops = dist_ops)
  
  # Stop the cluster
  stopCluster(cl)
  
  # Combine results and remove unnecessary columns
  df_tp_FL1_PSA <- do.call(rbind, l_PSA_input_FL1)
  df_tp_FL1_PSA[, -(1:7)]
}


# Make rank plot: function to visualize outcomes
# Create data sets for plotting: a new variable column for each line of treatment
make_rank_plot <- function(data, rank, outcome, savename, plottitle, scenario, top){
  
  sorted_data          <- data[order(rank),]       # step 1: sort by NHB
  sorted_data$rank_NHB <- seq(1:nrow(sorted_data)) # step 2: create rank vector
  top_sorted_data      <- sorted_data[1:top,]       # step 3: subset top 100
  table_data           <- sorted_data[1:10, c("Treatment_sequence",outcome)] # step 4: subset top 10
  
  sorted_data_rankCE   <- data[order(-data$"NHB"),]
  table_data_rankCE    <- sorted_data_rankCE[1:10, c("Treatment_sequence", "Cost", "QALYs", "NHB", "Time_in_remission")] 
  
  
  # Oude volgorde
  rank_sorted          <- data.frame("AZA"   = colSums(top_sorted_data[,2:6] == "AZA"), #MV EDIT 4 columns not 5
                                     "MTX"    = colSums(top_sorted_data[,2:6] == "MTX"),
                                     "IFX+AZA"    = colSums(top_sorted_data[,2:6] == "IFX+AZA"), 
                                     "IFX5"    = colSums(top_sorted_data[,2:6] == "IFX5"), 
                                     "ADA40"    = colSums(top_sorted_data[,2:6] == "ADA40"),
                                     "UST"    = colSums(top_sorted_data[,2:6] == "UST"), 
                                     "RIS"    = colSums(top_sorted_data[,2:6] == "RIS"),
                                     "UPA"    = colSums(top_sorted_data[,2:6] == "UPA"), 
                                     "VED"    = colSums(top_sorted_data[,2:6] == "VED")
                                      )
  
  rank_sorted_p      <- rank_sorted/ rowSums(rank_sorted) # df counts as proportion
  rank_sorted_p$line <- as.factor(c(paste("Line 1 \n", round(mean(top_sorted_data$Prop_in_line1, na.rm = T)), "% ind \n", 
                                          round(mean(top_sorted_data$Prop_mai_in_line1, na.rm = T)), "% mai \n", 
                                          round(mean(top_sorted_data$Time_in_line1, na.rm = T)),"yr mai",  sep = ""),
                                    paste("Line 2 \n", round(mean(top_sorted_data$Prop_in_line2, na.rm = T)), "% ind \n", 
                                          round(mean(top_sorted_data$Prop_mai_in_line2, na.rm = T)), "% mai \n", 
                                          round(mean(top_sorted_data$Time_in_line2, na.rm = T)),"yr mai",  sep = ""), 
                                    paste("Line 3 \n", round(mean(top_sorted_data$Prop_in_line3, na.rm = T)), "% ind \n", 
                                          round(mean(top_sorted_data$Prop_mai_in_line3, na.rm = T)), "% mai \n", 
                                          round(mean(top_sorted_data$Time_in_line3, na.rm = T)),"yr mai",  sep = ""), 
                                    paste("Line 4 \n", round(mean(top_sorted_data$Prop_in_line4, na.rm = T)), "% ind \n", 
                                          round(mean(top_sorted_data$Prop_mai_in_line4, na.rm = T)), "% mai \n", 
                                          round(mean(top_sorted_data$Time_in_line4, na.rm = T)),"yr mai",  sep = ""), 
                                    paste("Line 5 \n", round(mean(top_sorted_data$Prop_in_line5, na.rm = T)), "% ind \n", 
                                          round(mean(top_sorted_data$Prop_mai_in_line5, na.rm = T)), "% mai \n", 
                                          round(mean(top_sorted_data$Time_in_line5, na.rm = T)),"yr mai",  sep = "")))
  rank_sorted_p_long <- melt(rank_sorted_p)
  
  plot <- ggplot(data = rank_sorted_p_long, aes(x = line, y = value, fill = variable)) + 
    geom_bar(stat = "identity") +
    xlab("Treatment line") +
    #ylab(paste("Proportion top", top)) +
    ylab(paste("Proportion top 20%")) +
    scale_y_continuous(breaks = seq(0, 1, 0.1)) +
    scale_fill_manual(values=  c("AZA"          = "#5F6A6A" , 
                                 "MTX"          = "#85C1E9" ,
                                 "IFX.AZA"      = "#EBDEF0" ,
                                 "IFX5"          = "#C39BD3" , 
                                 "ADA40"        = "#76448A", 
                                 "UST"          = "#D4EFDF" , 
                                 "RIS"          = "#7DCEA0" ,
                                 "VED"          = "#EDBB99", 
                                 "UPA"          = "#ffba01"   
    ))+ 
    geom_text(data = subset(rank_sorted_p_long, value !=0), aes(label = paste0(round(value*100),"%")), 
              position = position_stack(vjust = 0.5), size = 2)+
    
    
    theme_bw() +
    ggtitle(paste(
      plottitle, 
      "\n Chart represents how often drug is in the top 20%", #"\n Scenario:",
      #"\n Chart represents how often drug is in the top", top, #"\n Scenario:",
      scenario)) +
    theme(legend.position = 'right',legend.title = element_blank())
  
  ggsave(here::here("visual output", 
                    paste(savename,  ".jpg", sep = "")), 
         width = 20, 
         height = 12, 
         units = c("cm"), 
         dpi = 300)
  
  write.table(table_data, file = here::here("visual output", 
                                            paste(savename, ".txt", sep = "")),
              sep  = ";",
              quote = FALSE, 
              row.names = F)  
  

  write.table(table_data_rankCE, file = here::here("visual output", 
                                                   paste(savename, "rankCE", ".txt", sep = "")),
              sep  = ";",
              quote = FALSE, 
              row.names = F)  
  
  output <- list(plot = plot, table_data = table_data, table_data_rankCE = table_data_rankCE)              
  
  return(output)
  
}


# Make rank plot: function to visualize outcomes
# Create data sets for plotting: a new variable column for each line of treatment
make_rank_plot_psa <- function(data, rank, outcome, savename, plottitle, scenario, top){
  
  sorted_data          <- data[order(rank),]       # step 1: sort by NHB
  sorted_data$rank_NHB <- seq(1:nrow(sorted_data)) # step 2: create rank vector
  top_sorted_data      <- sorted_data[1:top,]       # step 3: subset top 100
  table_data           <- sorted_data[1:10, c("Treatment_sequence",outcome)] # step 4: subset top 10
  
  sorted_data_rankCE   <- data[order(-data$"NHB"),]
  table_data_rankCE    <- sorted_data_rankCE[1:10, c("Treatment_sequence", "Cost", "QALYs", "NHB")] 
  
  
  # Oude volgorde
  rank_sorted          <- data.frame("AZA"   = colSums(top_sorted_data[,2:6] == "AZA"), #MV EDIT 4 columns not 5
                                     "MTX"    = colSums(top_sorted_data[,2:6] == "MTX"),
                                     "IFX+AZA"    = colSums(top_sorted_data[,2:6] == "IFX+AZA"), 
                                     "IFX5"    = colSums(top_sorted_data[,2:6] == "IFX5"), 
                                     "ADA40"    = colSums(top_sorted_data[,2:6] == "ADA40"),
                                     "UST"    = colSums(top_sorted_data[,2:6] == "UST"), 
                                     "RIS"    = colSums(top_sorted_data[,2:6] == "RIS"),
                                     "UPA"    = colSums(top_sorted_data[,2:6] == "UPA"), 
                                     "VED"    = colSums(top_sorted_data[,2:6] == "VED")
  )
 
  rank_sorted_p      <- rank_sorted/ rowSums(rank_sorted) # df counts as proportion
  rank_sorted_p$line <- as.factor(c(paste("Line 1"),
                                    paste("Line 2"), 
                                    paste("Line 3"), 
                                    paste("Line 4"), 
                                    paste("Line 5")))
    rank_sorted_p_long <- melt(rank_sorted_p)
  
  plot <- ggplot(data = rank_sorted_p_long, aes(x = line, y = value, fill = variable)) + 
    geom_bar(stat = "identity") +
    xlab("Treatment line") +
    #ylab(paste("Proportion top", top)) +
    ylab(paste("Proportion top 20%")) +
    scale_y_continuous(breaks = seq(0, 1, 0.1)) +
    scale_fill_manual(values=  c("AZA"          = "#5F6A6A" , 
                                 "MTX"          = "#85C1E9" ,
                                 "IFX.AZA"      = "#EBDEF0" ,
                                 "IFX5"          = "#C39BD3" , 
                                 "ADA40"        = "#76448A", 
                                 "UST"          = "#D4EFDF" , 
                                 "RIS"          = "#7DCEA0" ,
                                 "VED"          = "#EDBB99", 
                                 "UPA"          = "#ffba01"   
    ))+ 
    
    theme_bw() +
    ggtitle(paste(
      plottitle, 
      "\n Chart represents how often drug is in the top 20%", #"\n Scenario:",
      #"\n Chart represents how often drug is in the top", top, #"\n Scenario:",
      scenario)) +
    theme(legend.position = 'right',legend.title = element_blank())
  
  ggsave(here::here("visual output", 
                    paste(savename,  ".jpg", sep = "")), 
         width = 20, 
         height = 12, 
         units = c("cm"), 
         dpi = 300)
  
  write.table(table_data, file = here::here("visual output", 
                                            paste(savename, ".txt", sep = "")),
              sep  = ";",
              quote = FALSE, 
              row.names = F)  
  
  
  write.table(table_data_rankCE, file = here::here("visual output", 
                                                   paste(savename, "rankCE", ".txt", sep = "")),
              sep  = ";",
              quote = FALSE, 
              row.names = F)  
  
  output <- list(plot = plot, table_data = table_data, table_data_rankCE = table_data_rankCE)              
  
  return(output)
  
}
