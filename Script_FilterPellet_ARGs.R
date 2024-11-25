#### Import dataset ####

library(readxl)

Dataset_FIlterPellet <- read_excel("C:/Users/lenaf/Documents/R/Dataset_FIlterPellet.xlsx")
View(Dataset_FIlterPellet)


# Load necessary libraries
install.packages("tidyverse")
library(tidyverse)
install.packages("dplyr")
library(dplyr)
library(broom)
library(dplyr)
library(car)

# For analysis purposes, create a dataframe 'data'based on the dataset
data <- Dataset_FIlterPellet

#### mefA ####

# Selecting data for the 'mefA' gene
filter_mefA <- data %>% filter(`Concentration type` == 'filter', !is.na(mefA)) %>% pull(mefA)
pellet_mefA <- data %>% filter(`Concentration type` == 'pellet', !is.na(mefA)) %>% pull(mefA)

# Shapiro-Wilk Test for normality
shapiro_filter_mefA <- shapiro.test(filter_mefA)
shapiro_pellet_mefA <- shapiro.test(pellet_mefA)

# Printing Shapiro-Wilk results
print(shapiro_filter_mefA)
print(shapiro_pellet_mefA)

# Levene Test for equality of variances (if necessary)
if (shapiro_filter_mefA$p.value > 0.05 && shapiro_pellet_mefA$p.value > 0.05) {
  levene_mefA <- car::leveneTest(mefA ~ `Concentration type`, data = data)
  print(levene_mefA)
}

# T-Test or Mann-Whitney U Test (depending on previous results)
if (shapiro_filter_mefA$p.value > 0.05 && shapiro_pellet_mefA$p.value > 0.05 && levene_mefA$p.value > 0.05) {
  ttest_mefA <- t.test(mefA ~ `Concentration type`, data = data, var.equal = TRUE)
  print(ttest_mefA)
} else {
  mwtest_mefA <- wilcox.test(mefA ~ `Concentration type`, data = data)
  print(mwtest_mefA)
}

# Calculating and printing means
mean_filter_mefA <- mean(filter_mefA, na.rm = TRUE)
mean_pellet_mefA <- mean(pellet_mefA, na.rm = TRUE)
print(c(Mean_Filter = mean_filter_mefA, Mean_Pellet = mean_pellet_mefA))

#### mphE ####

# Selecting data for the 'mphE' gene
filter_mphE <- data %>% filter(`Concentration type` == 'filter', !is.na(mphE)) %>% pull(mphE)
pellet_mphE <- data %>% filter(`Concentration type` == 'pellet', !is.na(mphE)) %>% pull(mphE)

# Shapiro-Wilk Test for normality
shapiro_filter_mphE <- shapiro.test(filter_mphE)
shapiro_pellet_mphE <- shapiro.test(pellet_mphE)

# Printing Shapiro-Wilk results
print(shapiro_filter_mphE)
print(shapiro_pellet_mphE)

# Levene Test for equality of variances (if necessary)
if (shapiro_filter_mphE$p.value > 0.05 && shapiro_pellet_mphE$p.value > 0.05) {
  levene_mphE <- car::leveneTest(mphE ~ `Concentration type`, data = data)
  print(levene_mphE)
}

# T-Test or Mann-Whitney U Test (depending on previous results)
if (shapiro_filter_mphE$p.value > 0.05 && shapiro_pellet_mphE$p.value > 0.05 && levene_mphE$p.value > 0.05) {
  ttest_mphE <- t.test(mphE ~ `Concentration type`, data = data, var.equal = TRUE)
  print(ttest_mphE)
} else {
  mwtest_mphE <- wilcox.test(mphE ~ `Concentration type`, data = data)
  print(mwtest_mphE)
}

# Calculating and printing means
mean_filter_mphE <- mean(filter_mphE, na.rm = TRUE)
mean_pellet_mphE <- mean(pellet_mphE, na.rm = TRUE)
print(c(Mean_Filter = mean_filter_mphE, Mean_Pellet = mean_pellet_mphE))


#### SHV ####

# Selecting data for the 'blaSHV_1' gene
filter_blaSHV_1 <- data %>% filter(`Concentration type` == 'filter', !is.na(blaSHV_1)) %>% pull(blaSHV_1)
pellet_blaSHV_1 <- data %>% filter(`Concentration type` == 'pellet', !is.na(blaSHV_1)) %>% pull(blaSHV_1)

# Shapiro-Wilk Test for normality
shapiro_filter_blaSHV_1 <- shapiro.test(filter_blaSHV_1)
shapiro_pellet_blaSHV_1 <- shapiro.test(pellet_blaSHV_1)

# Printing Shapiro-Wilk results
print(shapiro_filter_blaSHV_1)
print(shapiro_pellet_blaSHV_1)

# Levene Test for equality of variances (if necessary)
levene_blaSHV_1 <- car::leveneTest(blaSHV_1 ~ `Concentration type`, data = data)
print(levene_blaSHV_1)

# T-Test 
ttest_blaSHV_1 <- t.test(blaSHV_1 ~ `Concentration type`, data = data, var.equal = TRUE)
print(ttest_blaSHV_1)

# Calculating and printing means
mean_filter_blaSHV_1 <- mean(filter_blaSHV_1, na.rm = TRUE)
mean_pellet_blaSHV_1 <- mean(pellet_blaSHV_1, na.rm = TRUE)
print(c(Mean_Filter = mean_filter_blaSHV_1, Mean_Pellet = mean_pellet_blaSHV_1))


#### CTX-M ####

# Selecting data for the 'blaCTXM' gene
filter_blaCTXM <- data %>% filter(`Concentration type` == 'filter', !is.na(blaCTXM)) %>% pull(blaCTXM)
pellet_blaCTXM <- data %>% filter(`Concentration type` == 'pellet', !is.na(blaCTXM)) %>% pull(blaCTXM)

# Shapiro-Wilk Test for normality
shapiro_filter_blaCTXM <- shapiro.test(filter_blaCTXM)
shapiro_pellet_blaCTXM <- shapiro.test(pellet_blaCTXM)

# Printing Shapiro-Wilk results
print(shapiro_filter_blaCTXM)
print(shapiro_pellet_blaCTXM)

# Levene Test for equality of variances (if necessary)
levene_blaCTXM <- car::leveneTest(blaCTXM ~ `Concentration type`, data = data)
print(levene_blaCTXM)

# T-Test 
ttest_blaCTXM <- t.test(blaCTXM ~ `Concentration type`, data = data, var.equal = TRUE)
print(ttest_blaCTXM)

# Calculating and printing means
mean_filter_blaCTXM <- mean(filter_blaCTXM, na.rm = TRUE)
mean_pellet_blaCTXM <- mean(pellet_blaCTXM, na.rm = TRUE)
print(c(Mean_Filter = mean_filter_blaCTXM, Mean_Pellet = mean_pellet_blaCTXM))


#### NDM ####

# Selecting data for the 'blaNDM' gene
filter_blaNDM <- data %>% filter(`Concentration type` == 'filter', !is.na(blaNDM)) %>% pull(blaNDM)
pellet_blaNDM <- data %>% filter(`Concentration type` == 'pellet', !is.na(blaNDM)) %>% pull(blaNDM)

# Shapiro-Wilk Test for normality
shapiro_filter_blaNDM <- shapiro.test(filter_blaNDM)
shapiro_pellet_blaNDM <- shapiro.test(pellet_blaNDM)

# Printing Shapiro-Wilk results
print(shapiro_filter_blaNDM)
print(shapiro_pellet_blaNDM)

# Levene Test for equality of variances (if necessary)
if (shapiro_filter_blaNDM$p.value > 0.05 && shapiro_pellet_blaNDM$p.value > 0.05) {
  levene_blaNDM <- car::leveneTest(blaNDM ~ `Concentration type`, data = data)
  print(levene_blaNDM)
}

# T-Test or Mann-Whitney U Test (depending on previous results)
if (shapiro_filter_blaNDM$p.value > 0.05 && shapiro_pellet_blaNDM$p.value > 0.05 && levene_blaNDM$p.value > 0.05) {
  ttest_blaNDM <- t.test(blaNDM ~ `Concentration type`, data = data, var.equal = TRUE)
  print(ttest_blaNDM)
} else {
  mwtest_blaNDM <- wilcox.test(blaNDM ~ `Concentration type`, data = data)
  print(mwtest_blaNDM)
}

# Calculating and printing means
mean_filter_blaNDM <- mean(filter_blaNDM, na.rm = TRUE)
mean_pellet_blaNDM <- mean(pellet_blaNDM, na.rm = TRUE)
print(c(Mean_Filter = mean_filter_blaNDM, Mean_Pellet = mean_pellet_blaNDM)) 


#### TEM ####

# Selecting data for the 'blaTEM' gene
filter_blaTEM<- data %>% filter(`Concentration type` == 'filter', !is.na(blaTEM)) %>% pull(blaTEM)
pellet_blaTEM <- data %>% filter(`Concentration type` == 'pellet', !is.na(blaTEM)) %>% pull(blaTEM)

# Shapiro-Wilk Test for normality
shapiro_filter_blaTEM <- shapiro.test(filter_blaTEM)
shapiro_pellet_blaTEM <- shapiro.test(pellet_blaTEM)

# Printing Shapiro-Wilk results
print(shapiro_filter_blaTEM)
print(shapiro_pellet_blaTEM)

# Levene Test for equality of variances (if necessary)
if (shapiro_filter_blaTEM$p.value > 0.05 && shapiro_pellet_blaTEM$p.value > 0.05) {
  levene_blaTEM <- car::leveneTest(blaTEM ~ `Concentration type`, data = data)
  print(levene_blaTEM)
}

# T-Test or Mann-Whitney U Test (depending on previous results)
if (shapiro_filter_blaTEM$p.value > 0.05 && shapiro_pellet_blaTEM$p.value > 0.05 && levene_blaTEM$p.value > 0.05) {
  ttest_blaTEM <- t.test(blaTEM ~ `Concentration type`, data = data, var.equal = TRUE)
  print(ttest_blaTEM)
} else {
  mwtest_blaTEM <- wilcox.test(blaTEM ~ `Concentration type`, data = data)
  print(mwtest_blaTEM)
}

# Calculating and printing means
mean_filter_blaTEM <- mean(filter_blaTEM, na.rm = TRUE)
mean_pellet_blaTEM <- mean(pellet_blaTEM, na.rm = TRUE)
print(c(Mean_Filter = mean_filter_blaTEM, Mean_Pellet = mean_pellet_blaTEM))


#### OXA-1/blaOXA-30 ####

# Selecting data for the 'blaOXA' gene
filter_blaOXA<- data %>% filter(`Concentration type` == 'filter', !is.na(blaOXA)) %>% pull(blaOXA)
pellet_blaOXA <- data %>% filter(`Concentration type` == 'pellet', !is.na(blaOXA)) %>% pull(blaOXA)

# Shapiro-Wilk Test for normality
shapiro_filter_blaOXA <- shapiro.test(filter_blaOXA)
shapiro_pellet_blaOXA <- shapiro.test(pellet_blaOXA)

# Printing Shapiro-Wilk results
print(shapiro_filter_blaOXA)
print(shapiro_pellet_blaOXA)

# Levene Test for equality of variances (if necessary)
if (shapiro_filter_blaOXA$p.value > 0.05 && shapiro_pellet_blaOXA$p.value > 0.05) {
  levene_blaOXA <- car::leveneTest(blaOXA ~ `Concentration type`, data = data)
  print(levene_blaOXA)
}

# T-Test or Mann-Whitney U Test (depending on previous results)
if (shapiro_filter_blaOXA$p.value > 0.05 && shapiro_pellet_blaOXA$p.value > 0.05 && levene_blaOXA$p.value > 0.05) {
  ttest_blaOXA <- t.test(blaOXA ~ `Concentration type`, data = data, var.equal = TRUE)
  print(ttest_blaOXA)
} else {
  mwtest_blaOXA <- wilcox.test(blaOXA ~ `Concentration type`, data = data)
  print(mwtest_blaOXA)
}

# Calculating and printing means
mean_filter_blaOXA <- mean(filter_blaOXA, na.rm = TRUE)
mean_pellet_blaOXA <- mean(pellet_blaOXA, na.rm = TRUE)
print(c(Mean_Filter = mean_filter_blaOXA, Mean_Pellet = mean_pellet_blaOXA))


#### KPC ####

# Selecting data for the 'blaKPC' gene
filter_blaKPC<- data %>% filter(`Concentration type` == 'filter', !is.na(blaKPC)) %>% pull(blaKPC)
pellet_blaKPC <- data %>% filter(`Concentration type` == 'pellet', !is.na(blaKPC)) %>% pull(blaKPC)

# Shapiro-Wilk Test for normality
shapiro_filter_blaKPC <- shapiro.test(filter_blaKPC)
shapiro_pellet_blaKPC <- shapiro.test(pellet_blaKPC)

# Printing Shapiro-Wilk results
print(shapiro_filter_blaKPC)
print(shapiro_pellet_blaKPC)

# Levene Test for equality of variances (if necessary)
if (shapiro_filter_blaKPC$p.value > 0.05 && shapiro_pellet_blaKPC$p.value > 0.05) {
  levene_blaKPC <- car::leveneTest(blaKPC ~ `Concentration type`, data = data)
  print(levene_blaKPC)
}

# T-Test or Mann-Whitney U Test (depending on previous results)
if (shapiro_filter_blaKPC$p.value > 0.05 && shapiro_pellet_blaKPC$p.value > 0.05 && levene_blaKPC$p.value > 0.05) {
  ttest_blaKPC <- t.test(blaKPC ~ `Concentration type`, data = data, var.equal = TRUE)
  print(ttest_blaKPC)
} else {
  mwtest_blaKPC <- wilcox.test(blaKPC ~ `Concentration type`, data = data)
  print(mwtest_blaKPC)
}

# Calculating and printing means
mean_filter_blaKPC <- mean(filter_blaKPC, na.rm = TRUE)
mean_pellet_blaKPC <- mean(pellet_blaKPC, na.rm = TRUE)
print(c(Mean_Filter = mean_filter_blaKPC, Mean_Pellet = mean_pellet_blaKPC))


#### IMP ####

# Selecting data for the 'blaIMP_1' gene
filter_blaIMP_1<- data %>% filter(`Concentration type` == 'filter', !is.na(blaIMP_1)) %>% pull(blaIMP_1)
pellet_blaIMP_1 <- data %>% filter(`Concentration type` == 'pellet', !is.na(blaIMP_1)) %>% pull(blaIMP_1)

# Shapiro-Wilk Test for normality
shapiro_filter_blaIMP_1 <- shapiro.test(filter_blaIMP_1)
shapiro_pellet_blaIMP_1 <- shapiro.test(pellet_blaIMP_1)

# Printing Shapiro-Wilk results
print(shapiro_filter_blaIMP_1)
print(shapiro_pellet_blaIMP_1)

# Levene Test for equality of variances (if necessary)
if (shapiro_filter_blaIMP_1$p.value > 0.05 && shapiro_pellet_blaIMP_1$p.value > 0.05) {
  levene_blaIMP_1 <- car::leveneTest(blaIMP_1 ~ `Concentration type`, data = data)
  print(levene_blaIMP_1)
}

# T-Test or Mann-Whitney U Test (depending on previous results)
if (shapiro_filter_blaIMP_1$p.value > 0.05 && shapiro_pellet_blaIMP_1$p.value > 0.05 && levene_blaIMP_1$p.value > 0.05) {
  ttest_blaIMP_1 <- t.test(blaIMP_1 ~ `Concentration type`, data = data, var.equal = TRUE)
  print(ttest_blaIMP_1)
} else {
  mwtest_blaIMP_1 <- wilcox.test(blaIMP_1 ~ `Concentration type`, data = data)
  print(mwtest_blaIMP_1)
}

# Calculating and printing means
mean_filter_blaIMP_1 <- mean(filter_blaIMP_1, na.rm = TRUE)
mean_pellet_blaIMP_1 <- mean(pellet_blaIMP_1, na.rm = TRUE)
print(c(Mean_Filter = mean_filter_blaIMP_1, Mean_Pellet = mean_pellet_blaIMP_1))


#### qnrB ####

# Selecting data for the 'qnrB' gene
filter_qnrB<- data %>% filter(`Concentration type` == 'filter', !is.na(qnrB)) %>% pull(qnrB)
pellet_qnrB <- data %>% filter(`Concentration type` == 'pellet', !is.na(qnrB)) %>% pull(qnrB)

# Shapiro-Wilk Test for normality
shapiro_filter_qnrB <- shapiro.test(filter_qnrB)
shapiro_pellet_qnrB <- shapiro.test(pellet_qnrB)

# Printing Shapiro-Wilk results
print(shapiro_filter_qnrB)
print(shapiro_pellet_qnrB)

# Levene Test for equality of variances (if necessary)
if (shapiro_filter_qnrB$p.value > 0.05 && shapiro_pellet_qnrB$p.value > 0.05) {
  levene_qnrB <- car::leveneTest(qnrB ~ `Concentration type`, data = data)
  print(levene_qnrB)
}

# T-Test or Mann-Whitney U Test (depending on previous results)
if (shapiro_filter_qnrB$p.value > 0.05 && shapiro_pellet_qnrB$p.value > 0.05 && levene_qnrB$p.value > 0.05) {
  ttest_qnrB <- t.test(qnrB ~ `Concentration type`, data = data, var.equal = TRUE)
  print(ttest_qnrB)
} else {
  mwtest_qnrB <- wilcox.test(qnrB ~ `Concentration type`, data = data)
  print(mwtest_qnrB)
}

# Calculating and printing means
mean_filter_qnrB <- mean(filter_qnrB, na.rm = TRUE)
mean_pellet_qnrB <- mean(pellet_qnrB, na.rm = TRUE)
print(c(Mean_Filter = mean_filter_qnrB, Mean_Pellet = mean_pellet_qnrB))

#### qnrA ####

# Selecting data for the 'qnrA' gene
filter_qnrA<- data %>% filter(`Concentration type` == 'filter', !is.na(qnrA)) %>% pull(qnrA)
pellet_qnrA <- data %>% filter(`Concentration type` == 'pellet', !is.na(qnrA)) %>% pull(qnrA)

# Shapiro-Wilk Test for normality
shapiro_filter_qnrA <- shapiro.test(filter_qnrA)
shapiro_pellet_qnrA <- shapiro.test(pellet_qnrA)

# Printing Shapiro-Wilk results
print(shapiro_filter_qnrA)
print(shapiro_pellet_qnrA)

# Levene Test for equality of variances (if necessary)
if (shapiro_filter_qnrA$p.value > 0.05 && shapiro_pellet_qnrA$p.value > 0.05) {
  levene_qnrA <- car::leveneTest(qnrA ~ `Concentration type`, data = data)
  print(levene_qnrA)
}

# T-Test or Mann-Whitney U Test (depending on previous results)
if (shapiro_filter_qnrA$p.value > 0.05 && shapiro_pellet_qnrA$p.value > 0.05 && levene_qnrA$p.value > 0.05) {
  ttest_qnrA <- t.test(qnrA ~ `Concentration type`, data = data, var.equal = TRUE)
  print(ttest_qnrA)
} else {
  mwtest_qnrA <- wilcox.test(qnrA ~ `Concentration type`, data = data)
  print(mwtest_qnrA)
}

# Calculating and printing means
mean_filter_qnrA <- mean(filter_qnrA, na.rm = TRUE)
mean_pellet_qnrA <- mean(pellet_qnrA, na.rm = TRUE)
print(c(Mean_Filter = mean_filter_qnrA, Mean_Pellet = mean_pellet_qnrA))