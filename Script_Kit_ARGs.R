#### Import dataset ####

library(readxl)

Dataset_Kit_ARGs <- read_excel("C:/Users/lenaf/Documents/R/Dataset_Kit_ARGs.xlsx")
View(Dataset_Kit_ARGs)


# Load necessary libraries
install.packages("tidyverse")
library(tidyverse)
install.packages("dplyr")
library(dplyr)
library(broom)
library(dplyr)
library(car)

# For analysis purposes, create a dataframe 'data'based on the dataset
data <- Dataset_Kit_ARGs

#### mefA ####

# Selecting data for the 'mefA' gene
PWP_mefA <- data %>% filter(`Kit` == 'PWP', !is.na(mefA)) %>% pull(mefA)
PL_mefA <- data %>% filter(`Kit` == 'PL', !is.na(mefA)) %>% pull(mefA)

# Calculating and printing means
mean_PWP_mefA <- mean(PWP_mefA, na.rm = TRUE)
mean_PL_mefA <- mean(PL_mefA, na.rm = TRUE)
print(c(Mean_PWP = mean_PWP_mefA, Mean_PL = mean_PL_mefA))

#### mphE ####

# Selecting data for the 'mphE' gene
PWP_mphE <- data %>% filter(`Kit` == 'PWP', !is.na(mphE)) %>% pull(mphE)
PL_mphE <- data %>% filter(`Kit` == 'PL', !is.na(mphE)) %>% pull(mphE)

# Calculating and printing means
mean_PWP_mphE <- mean(PWP_mphE, na.rm = TRUE)
mean_PL_mphE <- mean(PL_mphE, na.rm = TRUE)
print(c(Mean_PWP = mean_PWP_mphE, Mean_PL = mean_PL_mphE))


#### SHV ####

# Selecting data for the 'blaSHV_1' gene
PWP_blaSHV_1 <- data %>% filter(`Kit` == 'PWP', !is.na(blaSHV_1)) %>% pull(blaSHV_1)
PL_blaSHV_1 <- data %>% filter(`Kit` == 'PL', !is.na(blaSHV_1)) %>% pull(blaSHV_1)

# Calculating and printing means
mean_PWP_blaSHV_1 <- mean(PWP_blaSHV_1, na.rm = TRUE)
mean_PL_blaSHV_1 <- mean(PL_blaSHV_1, na.rm = TRUE)
print(c(Mean_PWP = mean_PWP_blaSHV_1, Mean_PL = mean_PL_blaSHV_1))


#### CTX-M ####

# Selecting data for the 'blaCTXM' gene
PWP_blaCTXM <- data %>% filter(`Kit` == 'PWP', !is.na(blaCTXM)) %>% pull(blaCTXM)
PL_blaCTXM <- data %>% filter(`Kit` == 'PL', !is.na(blaCTXM)) %>% pull(blaCTXM)

# Calculating and printing means
mean_PWP_blaCTXM <- mean(PWP_blaCTXM, na.rm = TRUE)
mean_PL_blaCTXM <- mean(PL_blaCTXM, na.rm = TRUE)
print(c(Mean_PWP = mean_PWP_blaCTXM, Mean_PL = mean_PL_blaCTXM))


#### NDM ####

# Selecting data for the 'blaNDM' gene
PWP_blaNDM <- data %>% filter(`Kit` == 'PWP', !is.na(blaNDM)) %>% pull(blaNDM)
PL_blaNDM <- data %>% filter(`Kit` == 'PL', !is.na(blaNDM)) %>% pull(blaNDM)

# Calculating and printing means
mean_PWP_blaNDM <- mean(PWP_blaNDM, na.rm = TRUE)
mean_PL_blaNDM <- mean(PL_blaNDM, na.rm = TRUE)
print(c(Mean_PWP = mean_PWP_blaNDM, Mean_PL = mean_PL_blaNDM)) 


#### TEM ####

# Selecting data for the 'blaTEM' gene
PWP_blaTEM<- data %>% filter(`Kit` == 'PWP', !is.na(blaTEM)) %>% pull(blaTEM)
PL_blaTEM <- data %>% filter(`Kit` == 'PL', !is.na(blaTEM)) %>% pull(blaTEM)

# Calculating and printing means
mean_PWP_blaTEM <- mean(PWP_blaTEM, na.rm = TRUE)
mean_PL_blaTEM <- mean(PL_blaTEM, na.rm = TRUE)
print(c(Mean_PWP = mean_PWP_blaTEM, Mean_PL = mean_PL_blaTEM))


#### OXA-1/blaOXA-30 ####

# Selecting data for the 'blaOXA' gene
PWP_blaOXA<- data %>% filter(`Kit` == 'PWP', !is.na(blaOXA)) %>% pull(blaOXA)
PL_blaOXA <- data %>% filter(`Kit` == 'PL', !is.na(blaOXA)) %>% pull(blaOXA)

# Calculating and printing means
mean_PWP_blaOXA <- mean(PWP_blaOXA, na.rm = TRUE)
mean_PL_blaOXA <- mean(PL_blaOXA, na.rm = TRUE)
print(c(Mean_PWP = mean_PWP_blaOXA, Mean_PL = mean_PL_blaOXA))


#### KPC ####

# Selecting data for the 'blaKPC' gene
PWP_blaKPC<- data %>% filter(`Kit` == 'PWP', !is.na(blaKPC)) %>% pull(blaKPC)
PL_blaKPC <- data %>% filter(`Kit` == 'PL', !is.na(blaKPC)) %>% pull(blaKPC)

# Calculating and printing means
mean_PWP_blaKPC <- mean(PWP_blaKPC, na.rm = TRUE)
mean_PL_blaKPC <- mean(PL_blaKPC, na.rm = TRUE)
print(c(Mean_PWP = mean_PWP_blaKPC, Mean_PL = mean_PL_blaKPC))


#### IMP ####

# Selecting data for the 'blaIMP_1' gene
PWP_blaIMP_1<- data %>% filter(`Kit` == 'PWP', !is.na(blaIMP_1)) %>% pull(blaIMP_1)
PL_blaIMP_1 <- data %>% filter(`Kit` == 'PL', !is.na(blaIMP_1)) %>% pull(blaIMP_1)

# Calculating and printing means
mean_PWP_blaIMP_1 <- mean(PWP_blaIMP_1, na.rm = TRUE)
mean_PL_blaIMP_1 <- mean(PL_blaIMP_1, na.rm = TRUE)
print(c(Mean_PWP = mean_PWP_blaIMP_1, Mean_PL = mean_PL_blaIMP_1))


#### qnrB ####

# Selecting data for the 'qnrB' gene
PWP_qnrB<- data %>% filter(`Kit` == 'PWP', !is.na(qnrB)) %>% pull(qnrB)
PL_qnrB <- data %>% filter(`Kit` == 'PL', !is.na(qnrB)) %>% pull(qnrB)

# Calculating and printing means
mean_PWP_qnrB <- mean(PWP_qnrB, na.rm = TRUE)
mean_PL_qnrB <- mean(PL_qnrB, na.rm = TRUE)
print(c(Mean_PWP = mean_PWP_qnrB, Mean_PL = mean_PL_qnrB))

#### qnrA ####

# Selecting data for the 'qnrA' gene
PWP_qnrA<- data %>% filter(`Kit` == 'PWP', !is.na(qnrA)) %>% pull(qnrA)
PL_qnrA <- data %>% filter(`Kit` == 'PL', !is.na(qnrA)) %>% pull(qnrA)

# Calculating and printing means
mean_PWP_qnrA <- mean(PWP_qnrA, na.rm = TRUE)
mean_PL_qnrA <- mean(PL_qnrA, na.rm = TRUE)
print(c(Mean_PWP = mean_PWP_qnrA, Mean_PL = mean_PL_qnrA))