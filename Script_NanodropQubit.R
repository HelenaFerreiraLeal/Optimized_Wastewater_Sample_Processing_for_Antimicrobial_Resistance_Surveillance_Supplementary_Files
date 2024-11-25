
#### Import Dataset ####

library(readxl)
NanodropQubitDataSet <- read_excel("NanodropQubitDataSet.xlsx", 
                                   col_types = c("text", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric", 
                                                 "numeric"))
View(NanodropQubitDataSet)

##### Load necessary library for better plotting ##### 
library(ggplot2)

result <- cor.test(NanodropQubitDataSet$`Log 10 Qubit`, NanodropQubitDataSet$`Log 10 Nanodrop`)

# Print the correlation coefficient and the p-value
print(paste("Pearson Correlation Coefficient: ", result$estimate))

###### Calculate the difference between Qubit and Nanodrop measurements #####
NanodropQubitDataSet$Difference <- NanodropQubitDataSet$`Log 10 Qubit` - NanodropQubitDataSet$`Log 10 Nanodrop`

#### Create a scatter plot with a regression line using ggplot2####
ggplot(NanodropQubitDataSet, aes(x=`Log 10 Qubit`, y=`Log 10 Nanodrop`)) +
  geom_point(aes(fill=Difference), size=3, stroke=0.5, shape=21) +  # Fill points based on difference
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0) +  # Gradient fill based on Difference
  geom_smooth(method="lm", se=FALSE, color="black", size=0.5) +  # Black regression line
  theme_minimal() +
  labs(title="Correlation between Qubit and Nanodrop",
       x="Log 10 Qubit",
       y="Log 10 Nanodrop") +
  theme(legend.position="right")  # Adjust legend position if needed


# correlation coefficient calculated between Nanodrop and 16s
cor_coef_nano <- cor(NanodropQubitDataSet$`Log 10 Nanodrop`, NanodropQubitDataSet$`Log 10 16S rRNA`, method = "pearson")


# Create the plot with black borders around points and fill color based on `Log 10 16s rRNA`
ggplot(NanodropQubitDataSet, aes(x=`Log 10 Nanodrop`, y=`Log 10 16S rRNA`)) +
  geom_point(aes(fill=`Log 10 16S rRNA`), size=3, stroke=0.5, shape=21) +  # Points with black border and colored fill
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=median(NanodropQubitDataSet$`Log 10 16S rRNA`), guide = 'colourbar') +  # Color scale with color guide
  geom_smooth(method="lm", se=FALSE, color="black", size=0.5) +  # Regression line
  theme_minimal() +
  labs(title="Correlation between Nanodrop and Log 10 16s rRNA",
       x="Log 10 Nanodrop",
       y="Log 10 16s rRNA") +
  annotate("text", x=Inf, y=Inf, label=paste("Pearson R:", round(cor_coef_nano, 2)), 
           hjust=1.1, vjust=1.1, size=5, color="black") +  # Correlation coefficient
  theme(legend.position="right")  # Legend position


#correlation coefficient calculated for Qubit and Log 10 16s rRNA
cor_coef_qubit <- cor(NanodropQubitDataSet$`Log 10 Qubit`, NanodropQubitDataSet$`Log 10 16S rRNA`, method = "pearson")

# Create the plot for Qubit with black borders around points and fill color based on `Log 10 16s rRNA`
ggplot(NanodropQubitDataSet, aes(x=`Log 10 Qubit`, y=`Log 10 16S rRNA`)) +
  geom_point(aes(fill=`Log 10 16S rRNA`), size=3, stroke=0.5, shape=21) +  # Points with black border and colored fill
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=median(NanodropQubitDataSet$`Log 10 16S rRNA`), guide = 'colourbar') +  # Color scale with color guide
  geom_smooth(method="lm", se=FALSE, color="black", size=0.5) +  # Regression line
  theme_minimal() +
  labs(title="Correlation between Qubit and Log 10 16s rRNA",
       x="Log 10 Qubit",
       y="Log 10 16s rRNA") +
  annotate("text", x=Inf, y=Inf, label=paste("Pearson R:", round(cor_coef_qubit, 2)), 
           hjust=1.1, vjust=1.1, size=5, color="black") +  # Correlation coefficient
  theme(legend.position="right")  # Legend position

# Correlation and P-value between Qubit and qPCR
result_qubit_qpcr <- cor.test(NanodropQubitDataSet$`Log 10 16S rRNA`, NanodropQubitDataSet$`Log 10 Qubit`, use="complete.obs")

# Print the correlation coefficient and the p-value for Qubit and qPCR
print(paste("Pearson Correlation Coefficient (Qubit vs qPCR): ", result_qubit_qpcr$estimate))
print(paste("P-value (Qubit vs qPCR): ", result_qubit_qpcr$p.value))

# Correlation and P-value between Nanodrop and qPCR
result_nanodrop_qpcr <- cor.test(NanodropQubitDataSet$`Log 10 16S rRNA`, NanodropQubitDataSet$`Log 10 Nanodrop`, use="complete.obs")

# Print the correlation coefficient and the p-value for Nanodrop and qPCR
print(paste("Pearson Correlation Coefficient (Nanodrop vs qPCR): ", result_nanodrop_qpcr$estimate))
print(paste("P-value (Nanodrop vs qPCR): ", result_nanodrop_qpcr$p.value))

