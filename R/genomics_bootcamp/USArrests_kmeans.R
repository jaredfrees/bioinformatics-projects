install.packages("tidyverse", dependencies = TRUE)

library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization

df <- USArrests

df <- na.omit(df)
head(df)

df <- scale(df)
head(df)


# default distance computed is the Euclidean
distance <- get_dist(df)
distance
# Heatmap
fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))


### k-means
k2 <- kmeans(df, centers = 2, nstart = 25)
k2

fviz_cluster(k2, data = df)

k3 <- kmeans(df, centers = 3, nstart = 25)
k4 <- kmeans(df, centers = 4, nstart = 25)
k5 <- kmeans(df, centers = 5, nstart = 25)

# plots to compare
p1 <- fviz_cluster(k2, geom = "point", data = df) + ggtitle("k = 2")
p2 <- fviz_cluster(k3, geom = "point",  data = df) + ggtitle("k = 3")
p3 <- fviz_cluster(k4, geom = "point",  data = df) + ggtitle("k = 4")
p4 <- fviz_cluster(k5, geom = "point",  data = df) + ggtitle("k = 5")

library(gridExtra)
grid.arrange(p1, p2, p3, p4, nrow = 2)

# Determine optimine number of k
set.seed(223456536)

# function to compute total within-cluster sum of square 
wss <- function(k) {
  kmeans(df, k, nstart = 10 )$tot.withinss
}

# Compute and plot wss for k = 1 to k = 15
k.values <- 1:15

# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)




fviz_nbclust(df, kmeans, method = "wss")
fviz_nbclust(df, kmeans, method = "silhouette")

gap_stat <- clusGap(df, FUN = kmeans, nstart = 25,
                    K.max = 10, B = 50)
fviz_gap_stat(gap_stat)



kmeans_clust<-eclust(df, k=2)



### PCA
pca = prcomp(df)
pca
pca.ind <- get_pca_ind(pca)

eig <- get_eigenvalue(pca)
fviz_eig(pca, col.var="blue")
fviz_pca_ind(pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("darkorchid4", "gold", "darkorange"),
             repel = TRUE
)