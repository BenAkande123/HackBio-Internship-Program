# HackBio Internship Program - February, 2025
# Team-arginine on the stage-3 task
# The group consists of four people, namely:
# Meltem (GitHub: @meltemktn)
# BenAkande (GitHub: @BenAkande123)
# Manav Vaish (GitHub: @manavvaish)
# Favour_Imoniye (GitHub: @Favour-Imoniye)
# GitHub link to the team-arginine repository for stage-3:

#Install packages
install.packages("DescTools",dependencies = c("Depends","Imports"),lib= .libPaths()[1])
install.packages("irlba")
install.packages("ClusterR")
install.packages ("viridis")

#Load the libraries
library (dplyr)
library (ggplot2)
library(DescTools)
library(irlba)
library(ClusterR)
library(viridis)

#Load the data set
Chem_descript <- read.csv("https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/drug_class_struct.txt", sep = "\t", header = TRUE)

#Assess the dimensions and datatypes of the data and the distribution of missing values in the data set
dim (Chem_descript)
colSums(is.na(Chem_descript))
str(Chem_descript)

#Remove non_numeric columns and columns with constant data or 0 variance for the principal component analysis (PCA)
Chemdescript_num <- Chem_descript %>% select (-c(ID, SMILES, target, ComponentCount, PosCount,  NegCount))

#Winsorize the data to minimize the impact of outliers during the scaling process
Chemdescript_winsor <- as.data.frame(lapply(Chemdescript_num, function(x) {
  Winsorize(x, val = quantile(x, probs = c(0.01, 0.99), na.rm = TRUE))})) #Calculates and caps all values at the 1st and 99th percentiles

#Check for the variance of each column
apply(Chemdescript_winsor, 2, var)

#Scale data to enhance analyses performance
Chemdescript_scaled <- scale(Chemdescript_winsor)
head(Chemdescript_scaled)

#Visualize the distribution of outliers in the data
boxplot(Chemdescript_scaled, las = 2, col = "lightgray", pch = 16, cex = 0.6,
        main = "Boxplot of Molecular Properties")

# Compute the PCA
Chem_pca <- prcomp(Chemdescript_scaled, center = TRUE, scale. = TRUE, rank. = 8)  # Keep only first 8 PCs

# Check explained variance
summary(Chem_pca)

# Select first N PCs explaining 90% variance of the overall data
cumulative_variance <- cumsum(Chem_pca$sdev^2) / sum(Chem_pca$sdev^2)
num_components <- which(cumulative_variance >= 0.90)[1]
Chempca_data <- Chem_pca$x[, 1:num_components]  # Keep relevant PCs

#Apply the elbow method to PCA reduced data
set.seed(123)

# Create a function to compute the total within-cluster sum of squares (WCSS)
wcss <- function(k) {
  kmeans(Chempca_data, centers = k, nstart = 1, iter.max=100)$tot.withinss
}

# Try different values of K (1 to 10)
k_values <- 1:10
wcss_values <- sapply(k_values, wcss)

# Plot Elbow Curve
plot(k_values, wcss_values, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of Clusters K", ylab = "Total Within-Cluster Sum of Squares",
     main = "Elbow Method on PCA-Reduced Data")


# Run MiniBatch K-Means
chem_matrix <- as.matrix(data.frame(Chempca_data))
set.seed(123)
k <-5 #Select a suitable amount of clusters
kmeans_result <- MiniBatchKmeans(chem_matrix, clusters = k, batch_size = 1000, num_init = 3, max_iters =100 )

# Extract cluster assignments
assign_clusters <- function(point, centroids) {
  return(which.min(colSums((t(centroids) - point)^2)))  # Find nearest centroid
}

clusters <- apply(chem_matrix, 1, assign_clusters, centroids = kmeans_result$centroids)


# Visualize the chemical space by docking score
clusters <-as.factor(clusters)
ggplot(chem_matrix, aes(x = PC1, y = PC2, color = Chemdescript_scaled[,"score"])) +
  geom_point(alpha = 0.7, size = 1) +
  scale_color_viridis(option = "plasma", direction = -1) +
  theme_minimal()+
labs(title = "Docking scores across the chemical space", color ="Docking Score")


#Visualize the clusters by docking scores
ggplot(data= data.frame(chem_matrix, clusters, score= Chemdescript_scaled[,"score"]),aes(x = PC1, y = PC2, color = score, shape=clusters)) +
  geom_point(alpha = 0.7) +
  scale_color_gradient(low="darkblue",high="lightblue") +
  facet_wrap(~ clusters) +
  theme_minimal() +
  ggtitle("Docking Scores Across Clusters")

#Examine the number of compounds with low docking score
# Filter the compounds in the lowest scoring cluster
low_score_cluster <- data.frame(chem_matrix, clusters, score = Chemdescript_scaled[,"score"])

low_score_cluster_compounds <- low_score_cluster %>%
  filter(clusters == 2)

Lowscorecompounds <- low_score_cluster_compounds %>% filter (score <=-2)
dim(Lowscorecompounds)

#Interpretation of results
# Cluster 2 had the highest representation of low docking scores across all clusters in the chemical space
#It contains 3461 compounds
# However,there was no distinguishable sub cluster with a preferentially low number of docking score


