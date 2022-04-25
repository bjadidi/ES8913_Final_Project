## Final Project for ES8913 ##

#####
# Libraries ---------------------------------------------------------------


## Libraries

install.packages("readr")
library(tidyverse)
library(ggplot2)



# Read Data ---------------------------------------------------------------

rm(list=ls()) 

setwd("./Data")

cwd <- getwd()

# Read Data and save them in a List

file_list <- list.files(path = cwd)
data_list <- list()

for (i in 1:(length(file_list))){
  
      curr_file <- str_replace_all(string=paste("./", file_list[i]),
                                   pattern=" ", repl="")
      
      curr_data <- read_csv(curr_file, show_col_types = FALSE,col_names = true)[1:12]
      
      colnames(curr_data) <- c("Type","Num","X","Y","Z","Vx","Vy","Vz","Fx","Fy","Fz","D")
      

      data_list[[i]] <- curr_data 

}

setwd("../")

force_data <- read_csv("./Force_Data.csv")
power_data <- read_csv("./Power_Data.csv")
rsd_data <- read_csv("./RSD_Data.csv")
gt_data <- read_csv("./GT_Data.csv")


# Calculation of RSD ------------------------------------------------------

## Domain and Num of Bins ##
GridN <- c(4,2,4)
low_0_RSD <- c(-0.28,-0.15,-0.56)
max_0_RSD <- c(0.28,0.17,0.0)
##
GridNx <- GridN[1]
GridNy <- GridN[2]
GridNz <- GridN[3]

lowx0 <- low_0_RSD[1]
lowy0 <- low_0_RSD[2]
lowz0 <- low_0_RSD[3]

maxx0 <- max_0_RSD[1]
maxy0 <- max_0_RSD[2]
maxz0 <- max_0_RSD[3]


binlengthx <- (maxx0 - lowx0)/GridNx
binlengthy <- (maxy0 - lowy0)/GridNy
binlengthz <- (maxz0 - lowz0)/GridNz


RSD <- c()
time_list <- c()
Num_bins <- GridNx * GridNy * GridNz

for (time in 0:(length(data_list)-1)){
  cat("Calculation of RSD for time : ", time, "\n")
  t <- time + 1
  
  
  sum_var <- 0
  sum_mean <- 0
  
  np1T <- nrow(filter(data_list[[t]], Type == 1))
  
  np3T <- nrow(filter(data_list[[t]], Type == 3))
  
  P <- np1T / np3T
  
  lowx <- lowx0
  for (i_Nx in 1:GridNx){
    lowxn <- lowx + binlengthx
    
    lowy <- lowy0
    for (i_Ny in 1:GridNy){
      lowyn <- lowy + binlengthy 
      
      lowz <- lowz0
      for (i_Nz in 1:GridNz){
        lowzn <- lowz + binlengthz 
        ####
        np1 <- nrow(filter(data_list[[t]], Type == 1, lowx <= X, X < lowxn,
                                                     lowy <= Y, Y < lowyn, 
                                                     lowz <= Z, Z < lowzn))
        
        np3 <- nrow(filter(data_list[[t]], Type == 3, lowx <= X, X < lowxn,
                                                     lowy <= Y, Y < lowyn, 
                                                     lowz <= Z, Z < lowzn))
        
        if (np1 > 5 & np3> 5){
          
          sum_var <- (np1 / (np1+np3) - P) ^ 2 + sum_var
          sum_mean <- np1 / (np1+np3) + sum_mean
          
        }
        else{
          Num_bins <- Num_bins - 1
        }
        ###
        lowz <- lowzn
      }
      lowy <- lowyn
      
    }
    lowx <- lowxn
  }
  
  
  var <- sum_var / Num_bins
  mean <- sum_mean / Num_bins
  svar <- sqrt(var)
  
  RSD[t] <- svar/mean * 100
  
  if (RSD[t] > 100){
    RSD[t] <- 100
  }
  
  time_list[t] <- time

  
  cat("RSD : ", RSD[t], "\n")
  
}

Result_list <- list(time_list, RSD)



Result_df <- as.data.frame(Result_list)
colnames(Result_df) <- c("Time","RSD")


ggplot(data = Result_df, aes(x = Time, y = RSD)) + 
  geom_line() + 
  geom_point() + 
  labs(x = "Time (s)",
       y = "RSD (%)")

ggsave("./Results/RSD.png",width = 20, height = 20, units = "cm")

colnames(rsd_data) <- c("Time","RSD")
ggplot(data = rsd_data, aes(x = Time, y = RSD)) + 
  geom_line() + 
  geom_point() + 
  labs(x = "Time(s)",
       y = "RSD(%)")

ggsave("./Results/RSD_New.png",width = 20, height = 20, units = "cm")






# Calculation of Granular Temp. -------------------------------------------
colnames(gt_data) <- c("Z","GT")
ggplot(data = gt_data, aes(x = Z, y = GT)) + 
  geom_line() + 
  geom_point() + 
  labs(x = "Z(bin number)",
       y = "GT(m2/s2)")

ggsave("./Results/GT.png",width = 20, height = 20, units = "cm")

## Domain and Num of Bins ##
GridN <- c(1,1,4)
low_0_GT <- c(0.117,-0.14,-0.56)
max_0_GT <- c(0.137,0.14,0.0)
##
GridNx <- GridN[1]
GridNy <- GridN[2]
GridNz <- GridN[3]

lowx0 <- low_0_RSD[1]
lowy0 <- low_0_RSD[2]
lowz0 <- low_0_RSD[3]

maxx0 <- max_0_RSD[1]
maxy0 <- max_0_RSD[2]
maxz0 <- max_0_RSD[3]


binlengthx <- (maxx0 - lowx0)/GridNx
binlengthy <- (maxy0 - lowy0)/GridNy
binlengthz <- (maxz0 - lowz0)/GridNz

Num_bins <- GridNx * GridNy * GridNz
RSD <- c()
time_list <- c()

for (time in 0:(length(data_list)-1)){
  cat("Calculation of GT for time : ", time, "\n")
  t <- time + 1
  
  sum_var <- 0
  sum_mean <- 0
  
  np1T <- nrow(filter(data_list[[t]], Type == 1))
  
  np3T <- nrow(filter(data_list[[t]], Type == 3))
  
  P <- np1T / np3T
  
  lowx <- lowx0
  for (i_Nx in 1:GridNx){
    lowxn <- lowx + binlengthx
    
    lowy <- lowy0
    for (i_Ny in 1:GridNy){
      lowyn <- lowy + binlengthy 
      
      lowz <- lowz0
      for (i_Nz in 1:GridNz){
        lowzn <- lowz + binlengthz 
        ####
        np1 <- nrow(filter(data_list[[t]], Type == 1, lowx <= X, X < lowxn,
                           lowy <= Y, Y < lowyn, 
                           lowz <= Z, Z < lowzn))
        
        np3 <- nrow(filter(data_list[[t]], Type == 3, lowx <= X, X < lowxn,
                           lowy <= Y, Y < lowyn, 
                           lowz <= Z, Z < lowzn))
        
        if (np1 > 0 ){
          
          sum_var <- (np1 / (np1+np3) - P) ^ 2 + sum_var
          sum_mean <- np1 / (np1+np3) + sum_mean
          
        }
        else{
          Num_bins <- Num_bins - 1
        }
        ###
        lowz <- lowzn
      }
      lowy <- lowyn
      
    }
    lowx <- lowxn
  }
  
  
  var <- sum_var / Num_bins
  mean <- sum_mean / Num_bins
  svar <- sqrt(var)
  
  RSD[t] <- svar/mean * 100
  
  if (RSD[t] > 100){
    RSD[t] <- 100
  }
  
  time_list[t] <- time
  
  
  cat("GT : ", RSD[t], "\n")
  
}

Result_list <- list(time_list, RSD)



Result_df <- as.data.frame(Result_list)
colnames(Result_df) <- c("Time","RSD")


ggplot(data = Result_df, aes(x = Time, y = RSD)) + 
  geom_line() + 
  geom_point() + 
  labs(x = "Time (s)",
       y = "RSD (%)")

ggsave("./Results/RSD.png",width = 20, height = 20, units = "cm")






# Calculation of Velocity Profile -----------------------------------------
GridN <- c(1,1,4)
low_0_V <- c(0.117,-0.14,-0.56)
max_0_V <- c(0.137,0.14,0.0)
##
GridNx <- GridN[1]
GridNy <- GridN[2]
GridNz <- GridN[3]

lowx0 <- low_0_RSD[1]
lowy0 <- low_0_RSD[2]
lowz0 <- low_0_RSD[3]

maxx0 <- max_0_RSD[1]
maxy0 <- max_0_RSD[2]
maxz0 <- max_0_RSD[3]


binlengthx <- (maxx0 - lowx0)/GridNx
binlengthy <- (maxy0 - lowy0)/GridNy
binlengthz <- (maxz0 - lowz0)/GridNz


df <- data_list[[20]]
df$V <- (df$Vx^2 + df$Vx^2 + df$Vx^2)^0.5

slice_df <- filter(df, -0.01 <= Z, Z < 0)
slice_df <- filter(df, Z < min(Z) + 0.001, Z >= min(Z) - 0.001)

ggplot(slice_df, aes(x = X, y = Y)) +
  geom_segment(aes(xend = X , yend = Y, colour = V),
               arrow = arrow(length = unit(0.1, "cm")),size = 0.5)+
  scale_colour_continuous(low = "grey80", high = "darkred")

ggsave("./Results/Vel_Prof.png",width = 20, height = 20, units = "cm")


# Calculation of Force acting on particles --------------------------------
colnames(force_data) <- c("Z","F")
ggplot(data = force_data, aes(x = Z, y = F)) + 
  geom_line() + 
  geom_point() + 
  labs(x = "Z(bin number)",
       y = "Force(N)")

ggsave("./Results/Force.png",width = 20, height = 20, units = "cm")

# Calculation of Force and Torque acting on geometry ----------------------


# Calculation of system's power consumption --------------------------------
power_data$F <- (power_data$Fx^2 + power_data$Fy^2 + power_data$Fz^2)^ 0.5 
power_data$M <- (power_data$Mx^2 + power_data$My^2 + power_data$Mz^2)^ 0.5 
power_data$Time <- power_data$t/(1040*20)

ggplot(data = power_data, aes(x = Time, y = F)) + 
  geom_line() + 
  geom_point() + 
  labs(x = "Time(s)",
       y = "Force(N)")

ggsave("./Results/Force_Impellers.png",width = 20, height = 20, units = "cm")

ggplot(data = power_data, aes(x = Time, y = M)) + 
  geom_line() + 
  geom_point() + 
  labs(x = "Time(s)",
       y = "Torque(N.m)")

ggsave("./Results/Torque_Impellers.png",width = 20, height = 20, units = "cm")


# Calculation of diffusion coef. and Peclet number ----------------------
#####