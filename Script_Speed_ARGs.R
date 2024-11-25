#### Import dataset ####

library(readxl)

Dataset_Speed_ARGs <- read_excel("C:/Users/lenaf/Documents/R/Dataset_Speed_ARGs.xlsx")
View(Dataset_Speed_ARGs)

# Load necessary libraries
install.packages("tidyverse")
library(tidyverse)
install.packages("dplyr")
library(dplyr)
library(broom)
library(dplyr)
library(car)

# For analysis purposes, create a dataframe 'data'based on the dataset
data <- Dataset_Speed_ARGs

#### mefA ####

# Selecting data for the 'mefA' gene
Sp3400_mefA <- data %>% filter(`Centrifugation Speed` == '3400', !is.na(mefA)) %>% pull(mefA)
Sp15000_mefA <- data %>% filter(`Centrifugation Speed` == '15000', !is.na(mefA)) %>% pull(mefA)

# Calculating and printing means
mean_3400_mefA <- mean(Sp3400_mefA, na.rm = TRUE)
mean_15000_mefA <- mean(Sp15000_mefA, na.rm = TRUE)
print(c(Mean_3400 = mean_3400_mefA, Mean_15000 = mean_15000_mefA))

#### mphE ####

# Selecting data for the 'mphE' gene
Sp3400_mphE <- data %>% filter(`Centrifugation Speed` == '3400', !is.na(mphE)) %>% pull(mphE)
Sp15000_mphE <- data %>% filter(`Centrifugation Speed` == '15000', !is.na(mphE)) %>% pull(mphE)

# Calculating and printing means
mean_3400_mphE <- mean(Sp3400_mphE, na.rm = TRUE)
mean_15000_mphE <- mean(Sp15000_mphE, na.rm = TRUE)
print(c(Mean_3400 = mean_3400_mphE, Mean_15000 = mean_15000_mphE))


#### SHV ####

# Selecting data for the 'blaSHV_1' gene
Sp3400_blaSHV_1 <- data %>% filter(`Centrifugation Speed` == '3400', !is.na(blaSHV_1)) %>% pull(blaSHV_1)
Sp15000_blaSHV_1 <- data %>% filter(`Centrifugation Speed` == '15000', !is.na(blaSHV_1)) %>% pull(blaSHV_1)

# Calculating and printing means
mean_3400_blaSHV_1 <- mean(Sp3400_blaSHV_1, na.rm = TRUE)
mean_15000_blaSHV_1 <- mean(Sp15000_blaSHV_1, na.rm = TRUE)
print(c(Mean_3400 = mean_3400_blaSHV_1, Mean_15000 = mean_15000_blaSHV_1))


#### CTX-M ####

# Selecting data for the 'blaCTXM' gene
Sp3400_blaCTXM <- data %>% filter(`Centrifugation Speed` == '3400', !is.na(blaCTXM)) %>% pull(blaCTXM)
Sp15000_blaCTXM <- data %>% filter(`Centrifugation Speed` == '15000', !is.na(blaCTXM)) %>% pull(blaCTXM)

# Calculating and printing means
mean_3400_blaCTXM <- mean(Sp3400_blaCTXM, na.rm = TRUE)
mean_15000_blaCTXM <- mean(Sp15000_blaCTXM, na.rm = TRUE)
print(c(Mean_3400 = mean_3400_blaCTXM, Mean_15000 = mean_15000_blaCTXM))


#### NDM ####

# Selecting data for the 'blaNDM' gene
Sp3400_blaNDM <- data %>% filter(`Centrifugation Speed` == '3400', !is.na(blaNDM)) %>% pull(blaNDM)
Sp15000_blaNDM <- data %>% filter(`Centrifugation Speed` == '15000', !is.na(blaNDM)) %>% pull(blaNDM)

# Calculating and printing means
mean_3400_blaNDM <- mean(Sp3400_blaNDM, na.rm = TRUE)
mean_15000_blaNDM <- mean(Sp15000_blaNDM, na.rm = TRUE)
print(c(Mean_3400 = mean_3400_blaNDM, Mean_15000 = mean_15000_blaNDM)) 


#### TEM ####

# Selecting data for the 'blaTEM' gene
Sp3400_blaTEM<- data %>% filter(`Centrifugation Speed` == '3400', !is.na(blaTEM)) %>% pull(blaTEM)
Sp15000_blaTEM <- data %>% filter(`Centrifugation Speed` == '15000', !is.na(blaTEM)) %>% pull(blaTEM)

# Calculating and printing means
mean_3400_blaTEM <- mean(Sp3400_blaTEM, na.rm = TRUE)
mean_15000_blaTEM <- mean(Sp15000_blaTEM, na.rm = TRUE)
print(c(Mean_3400 = mean_3400_blaTEM, Mean_15000 = mean_15000_blaTEM))


#### OXA-1/blaOXA-30 ####

# Selecting data for the 'blaOXA' gene
Sp3400_blaOXA<- data %>% filter(`Centrifugation Speed` == '3400', !is.na(blaOXA)) %>% pull(blaOXA)
Sp15000_blaOXA <- data %>% filter(`Centrifugation Speed` == '15000', !is.na(blaOXA)) %>% pull(blaOXA)

# Calculating and printing means
mean_3400_blaOXA <- mean(Sp3400_blaOXA, na.rm = TRUE)
mean_15000_blaOXA <- mean(Sp15000_blaOXA, na.rm = TRUE)
print(c(Mean_3400 = mean_3400_blaOXA, Mean_15000 = mean_15000_blaOXA))


#### KPC ####

# Selecting data for the 'blaKPC' gene
Sp3400_blaKPC<- data %>% filter(`Centrifugation Speed` == '3400', !is.na(blaKPC)) %>% pull(blaKPC)
Sp15000_blaKPC <- data %>% filter(`Centrifugation Speed` == '15000', !is.na(blaKPC)) %>% pull(blaKPC)

# Calculating and printing means
mean_3400_blaKPC <- mean(Sp3400_blaKPC, na.rm = TRUE)
mean_15000_blaKPC <- mean(Sp15000_blaKPC, na.rm = TRUE)
print(c(Mean_3400 = mean_3400_blaKPC, Mean_15000 = mean_15000_blaKPC))


#### IMP ####

# Selecting data for the 'blaIMP_1' gene
Sp3400_blaIMP_1<- data %>% filter(`Centrifugation Speed` == '3400', !is.na(blaIMP_1)) %>% pull(blaIMP_1)
Sp15000_blaIMP_1 <- data %>% filter(`Centrifugation Speed` == '15000', !is.na(blaIMP_1)) %>% pull(blaIMP_1)

# Calculating and printing means
mean_3400_blaIMP_1 <- mean(Sp3400_blaIMP_1, na.rm = TRUE)
mean_15000_blaIMP_1 <- mean(Sp15000_blaIMP_1, na.rm = TRUE)
print(c(Mean_3400 = mean_3400_blaIMP_1, Mean_15000 = mean_15000_blaIMP_1))


#### qnrB ####

# Selecting data for the 'qnrB' gene
Sp3400_qnrB<- data %>% filter(`Centrifugation Speed` == '3400', !is.na(qnrB)) %>% pull(qnrB)
Sp15000_qnrB <- data %>% filter(`Centrifugation Speed` == '15000', !is.na(qnrB)) %>% pull(qnrB)

# Calculating and printing means
mean_3400_qnrB <- mean(Sp3400_qnrB, na.rm = TRUE)
mean_15000_qnrB <- mean(Sp15000_qnrB, na.rm = TRUE)
print(c(Mean_3400 = mean_3400_qnrB, Mean_15000 = mean_15000_qnrB))

#### qnrA ####

# Selecting data for the 'qnrA' gene
Sp3400_qnrA<- data %>% filter(`Centrifugation Speed` == '3400', !is.na(qnrA)) %>% pull(qnrA)
Sp15000_qnrA <- data %>% filter(`Centrifugation Speed` == '15000', !is.na(qnrA)) %>% pull(qnrA)

# Calculating and printing means
mean_3400_qnrA <- mean(Sp3400_qnrA, na.rm = TRUE)
mean_15000_qnrA <- mean(Sp15000_qnrA, na.rm = TRUE)
print(c(Mean_3400 = mean_3400_qnrA, Mean_15000 = mean_15000_qnrA))