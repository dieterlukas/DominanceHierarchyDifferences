
# Load the necessary libraries with the matrices and to calculate the Bayesian Steepness measure. To install these, see https://github.com/gobbios/EloSteepness.data

library(dplyr)
library(EloSteepness.data)
library(EloSteepness)
library(EloRating)

# Pick only the matrices for chimpanzees
chimpanzees<-empirical_raw_data[names(empirical_raw_data) %in% c("watts1998_dom_1","newton-fisher2004_dom_1","newton-fisher2004_dom_2","muller2004_dom_1","muller2004_dom_2","mitani1993_dom_1","funkhouser2018_dom_1","wittig2003_dom_1")]


# Create two empty vector to store the output for comparison. For each matrix in the subset, we want to know the steepness based on the original data, and the steepness calculated for a matrix with the reduced sampling.
originalsteepness<-NA
reducedsteepness<-NA

# For each of the 8 matrices in the chimpanzee subset, first calculate the original steepness. Next, reduce the sampling by reducing the number of observed interactions by the half the average number of interactions. In cases where this leads to a negative number, assume that no interactions have been observed. Then calculate the steepness for this matrix with the reduced number of interactions.
for (i in 1:length(names(chimpanzees))){
  currenthierarchy<-chimpanzees[[i]]
  currentoriginalsteepness <- elo_steepness_from_matrix(mat = currenthierarchy, 
                                                 n_rand = case_when(sum(currenthierarchy)>0 ~ 20,sum(currenthierarchy)>500 ~ 5,sum(currenthierarchy)<100 ~ 50  ), 
                                                 refresh = 0, 
                                                 cores = 2, 
                                                 iter = 1000, 
                                                 seed = 1)
  originalsteepness[i]<-median(currentoriginalsteepness$steepness)
  print(c("completed original",i))
  reducedhierarchy<-round(currenthierarchy/2,0)
  reducedhierarchy[reducedhierarchy<0]<-0
  currentreducedsteepness<- elo_steepness_from_matrix(mat = reducedhierarchy, 
                                               n_rand = case_when(sum(reducedhierarchy)>0 ~ 20,sum(reducedhierarchy)>500 ~ 5,sum(reducedhierarchy)<100 ~ 50  ), 
                                               refresh = 0, 
                                               cores = 2, 
                                               iter = 1000, 
                                               seed = 1)
  reducedsteepness[i]<-median(currentreducedsteepness$steepness)
  print(c("completed reduced",i))
}

plot(reducedsteepness~originalsteepness,xlim=c(0,1),ylim=c(0,1))
abline(a=0,b=1)

chimpestimates<-as.data.frame(matrix(ncol=2,nrow=8))
colnames(chimpestimates)<-c("originalsteepness","reducedsteepness")
chimpestimates$originalsteepness<-originalsteepness
chimpestimates$reducedsteepness<-reducedsteepness
rownames(chimpestimates)<-c("watts1998_dom_1","newton-fisher2004_dom_1","newton-fisher2004_dom_2","muller2004_dom_1","muller2004_dom_2","mitani1993_dom_1","funkhouser2018_dom_1","wittig2003_dom_1")
