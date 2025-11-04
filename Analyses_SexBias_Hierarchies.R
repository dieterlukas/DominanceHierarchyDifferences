# This is code is for the analyses of the article by Spicher, Huchard, Lukas 
# The aim is to identify whether the hierarchies of males and females differ among primates

# To address this aim, we will answer the following questions:

# A) Are our measures robust?  
#   1) Are hierarchy steepness and linearity sensitive to sparseness and nb of interactions? 
#   
#   2) Are hierarchy steepness and linearity sensitive to the nb of individuals included in the hierarchy?

# B) Are there dominance types?  
#   3) Are hierarchy steepness and linearity correlated across species (w. or without accounting for phylogenetic relatedness)?
#   4) Are hierarchies steeper and more linear when they are based on signals rather than aggression?

# C) Are there sex-differences in the characteristics of hierarchies?
#   5) Is hierarchy steepness in males different from that in females (w. or without accounting for phylogenetic relatedness)?
#   
#   6) If steepness differs in males and females, can it be linked to the fact that female hierarchies often include more individuals?
#   
#   7) Is hierarchy linearity in males different from that in females (w. or without accounting for phylogenetic relatedness)?
#   
#   8) Is hierarchy linearity linked to hierarchy steepness in the same way in males and in females?
#   
#   9) Are hierarchies steeper for the sex that wins more fights - that is, is the proportion of intersexual fights that females win negatively related to the hierarchy steepness in males, and positively to the hierarchy steepness in females?
#   
#   10) Are hierarchies for females more likely to be based on signals whereas those in males more likely to be based on aggression?


# We first load the required packages
library(rethinking)
library(dplyr)
library(ape)
library(geiger)
library(phytools)

# Load the data
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
data <- read.csv("Comparative_primate_steepness_values.csv")


# Load the phylogeny and prepare it

specieslist<-as.data.frame(matrix(unique(data$species),ncol=1,nrow=length(unique(data$species))))
colnames(specieslist)<-"species"
rownames(specieslist)<-specieslist$species

# Load the phylogeny
phylogeny <- read.nexus("Upham2019MammalPhylogeny.nex")

# check whether the species in our data are the same as the species in the phylogeny
speciesmatching<-name.check(phylogeny,specieslist)
# check which species from our data are not in the phylogeny
speciesmatching$data_not_tree

specieslist<-as.data.frame(matrix(unique(data$species),ncol=1,nrow=length(unique(data$species))))
colnames(specieslist)<-"species"
rownames(specieslist)<-specieslist$species

# and reduce the phylogeny to only include the species for which we have data
speciesmatching<-name.check(phylogeny,specieslist)
mtree<-drop.tip(phylogeny,speciesmatching$tree_not_data)








# Analyses for section A - robustness

# 1a) steepness sensitive to sparseness (check the sparseness for outliers first)
dat_list_steepness_sparseness <- list(
  steepness = as.numeric(data$steepness),  
  sparseness = as.numeric(data$sparseness)
  )
# We need to provide priors, our expectation of what these values might be. For the means, we could expect that they are somewhere around 0.5
# For the variance, we expect this to be larger than than zero (so we use the dexp function) and larger than one
# For the relationship with sparseness, we don't expect any effect a priori
m_steepness_sparseness <- ulam(
  alist(
    steepness ~ dbeta2(mean,variance),
    logit(mean) <-a + b*sparseness,
    a ~dnorm(0.5,1),   
    b ~dnorm(0,1),
    variance ~ dexp(10)
  ) , data=dat_list_steepness_sparseness , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)
# We check the results - we are interested in the effect measured in the factor b
precis(m_steepness_sparseness)
plot(data$steepness~data$sparseness)
# We want this more like a meta-analytical model: sparseness describes precision, with sparser data the values go closer to a mean rather than the actual value. 


# Study duration influences sparseness, but has not separate influence on the measurements
# We need to remove missing values
durationdata<-data[is.na(data$studyperiod_month)==F,]
dat_list_steepness_sparseness_duration <- list(
  steepness = as.numeric(durationdata$steepness),  
  sparseness = as.numeric(durationdata$sparseness),
  studyduration = as.numeric(durationdata$studyperiod_month)
)
# We need to provide priors, our expectation of what these values might be. For the means, we could expect that they are somewhere around 0.5
# For the variance, we expect this to be larger than than zero (so we use the dexp function) and larger than one
# For the relationship with sparseness, we don't expect any effect a priori
m_steepness_sparseness_duration <- ulam(
  alist(
    steepness ~ dbeta2(mean,variance),
    logit(mean) <-a + b*sparseness + c*studyduration,
    a ~dnorm(0.5,1),   
    b ~dnorm(0,1),
    c ~dnorm(0,1),
    variance ~ dexp(10)
  ) , data=dat_list_steepness_sparseness_duration , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)
# We check the results - we are interested in the effect measured in the factor b and c
precis(m_steepness_sparseness_duration)



# 2a) steepness sensitive to number of interactions
dat_list_steepness_numberofineractions <- list(
  steepness = as.numeric(data$steepness),  
  numberofineractions = standardize(data$numberofineractions)
)
# We need to provide priors, our expectation of what these values might be. For the means, we could expect that they are somewhere around 0.5
# For the variance, we expect this to be larger than than zero (so we use the dexp function) and larger than one
# For the relationship with numberofineractions, we don't expect any effect a priori
m_steepness_numberofineractions <- ulam(
  alist(
    steepness ~ dbeta2(mean,variance),
    logit(mean) <-a + b*numberofineractions,
    a ~dnorm(0.5,1),   
    b ~dnorm(0,1),
    variance ~ dexp(10)
  ) , data=dat_list_steepness_numberofineractions , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)
# We check the results - we are interested in the effect measured in the factor b
precis(m_steepness_numberofineractions)
plot(data$steepness~log(data$numberofineractions))
# Here, the relationship is linear, so we might want to include this as a linear predictor


# 1b) linearity sensitive to sparseness
dat_list_h_index_sparseness <- list(
  h_index = as.numeric(data$h_index)-0.0001,  
  sparseness = standardize(data$sparseness)
)
# We need to provide priors, our expectation of what these values might be. For the means, we could expect that they are somewhere around 0.5
# For the variance, we expect this to be larger than than zero (so we use the dexp function) and larger than one
# For the relationship with sparseness, we don't expect any effect a priori
m_h_index_sparseness <- ulam(
  alist(
    h_index ~ dbeta2(mean,variance),
    logit(mean) <-a + b*sparseness,
    a ~dnorm(0.5,1),   
    b ~dnorm(0,1),
    variance ~ dexp(10)
  ) , data=dat_list_h_index_sparseness , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)
# We check the results - we are interested in the effect measured in the factor b
precis(m_h_index_sparseness)
plot(data$h_index~data$sparseness)
# Here, the relationship is linear, so we want to include this as a linear predictor


# 2b) linearity sensitive to number of interactions
dat_list_h_index_numberofineractions <- list(
  h_index = as.numeric(data$h_index)-0.0001,  
  numberofineractions = standardize(data$numberofineractions)
)
# We need to provide priors, our expectation of what these values might be. For the means, we could expect that they are somewhere around 0.5
# For the variance, we expect this to be larger than than zero (so we use the dexp function) and larger than one
# For the relationship with numberofineractions, we don't expect any effect a priori
m_h_index_numberofineractions <- ulam(
  alist(
    h_index ~ dbeta2(mean,variance),
    logit(mean) <-a + b*numberofineractions,
    a ~dnorm(0.5,1),   
    b ~dnorm(0,1),
    variance ~ dexp(10)
  ) , data=dat_list_h_index_numberofineractions , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)
# We check the results - we are interested in the effect measured in the factor b
precis(m_h_index_numberofineractions)
plot(data$h_index~log(data$numberofineractions))
# Here is no effect, so we do not need to account for it


# 2c) linearity different in captivity than in the wild
dat_list_h_index_captivity <- list(
  h_index = as.numeric(data$h_index)-0.0001,  
  captive = as.numeric(as.factor(data$group)),
  numberofindividuals= standardize(data$numberofindividuals)
)
# We need to provide priors, our expectation of what these values might be. For the means, we could expect that they are somewhere around 0.5
# For the variance, we expect this to be larger than than zero (so we use the dexp function) and larger than one
# For the relationship with numberofineractions, we don't expect any effect a priori
m_h_index_captivity <- ulam(
  alist(
    h_index ~ dbeta2(mean,variance),
    logit(mean) <-a + b[captive]+c*numberofindividuals,
    a ~dnorm(0.5,1),   
    b[captive] ~dnorm(0,1),
    c ~dnorm(0,1),
    variance ~ dexp(10)
  ) , data=dat_list_h_index_captivity , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)
# We check the results - we are interested in the effect measured in the factor b
posterior_captivity<-extract.samples(m_h_index_captivity)
contrast_captivity<-inv_logit(posterior_captivity$b[,2])-inv_logit(posterior_captivity$b[,1])
precis(contrast_captivity)





### Analyses for section B - dominance types

#   3) Are hierarchy steepness and linearity correlated across species  / accounting for phylogenetic relatedness and sparseness (which affects both)

mdata_phylogeny_steepness_linearity <- list(
  steepness=data$steepness,
  linearity=data$h_index,
  sparseness=data$sparseness,
  species=as.integer(as.factor(data$species)),
  N_spp=length(unique(data$species))
)

Dmat<-cophenetic(mtree)
mdata_phylogeny_steepness_linearity$Dmat<-Dmat/max(Dmat)
colnames(mdata_phylogeny_steepness_linearity$Dmat)<-as.integer(as.factor(colnames(mdata_phylogeny_steepness_linearity$Dmat)))
rownames(mdata_phylogeny_steepness_linearity$Dmat)<-as.integer(as.factor(rownames(mdata_phylogeny_steepness_linearity$Dmat)))

m_steepness_steepness_linearity <- ulam(
  alist(
    steepness ~ dbeta2(mean,variance),
    logit(mean) <-a+b[species]+c*sparseness+d*linearity,
    a ~dnorm(0,1),
    c~dnorm(0,1),
    d~dnorm(0,1),
    vector[N_spp]:b~multi_normal(0,SIGMA),
    matrix[N_spp,N_spp]:SIGMA <- cov_GPL2( Dmat , etasq, rhosq , 0.01 ),
    etasq ~ half_normal(1,0.25),
    rhosq ~ half_normal(3,0.25),
    variance~dexp(10)
  ) , data=mdata_phylogeny_steepness_linearity , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)

precis(m_steepness_steepness_linearity)
plot(data$steepness~data$h_index)
# Linearity and steepness are positively correlated

#   4) Are hierarchies steeper and more linear when they are based on signals rather than aggression?
data[data$typeofbehaviour=="unc",]$typeofbehaviour<-"AD"
data[data$typeofbehaviour=="flee",]$typeofbehaviour<-"A"
boxplot(data$steepness~data$typeofbehaviour)
boxplot(data$h_index~data$typeofbehaviour)


dat_list_steepness_typeofbehaviour <- list(
  steepness = as.numeric(data$steepness),  
  typeofbehaviour = as.integer(as.factor(data$typeofbehaviour))
)
m_steepness_typeofbehaviour <- ulam(
  alist(
    steepness ~ dbeta2(mean,variance),
    logit(mean) <-a + b[typeofbehaviour],
    a ~dnorm(0.5,1),   
    b[typeofbehaviour] ~dnorm(0,1),
    variance ~ dexp(10)
  ) , data=dat_list_steepness_typeofbehaviour , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)
# We check the results - we are interested in the effect measured in the factor b
precis(m_steepness_typeofbehaviour)
post_steepness_typeofbehaviour<-extract.samples(m_steepness_typeofbehaviour)
# Comparing D versus A
contrast_D_A<-post_steepness_typeofbehaviour$b[,3]-post_steepness_typeofbehaviour$b[,1]
precis(contrast_D_A)
# No real difference - and AD is actually not in the middle, but the lowest.



### Analyses for section C - sex-differences 
#   5) Is hierarchy steepness in males different from that in females (w. or without accounting for phylogenetic relatedness)?

# descriptive statistics
length(unique(data$species)) # 38
table(data$sex) # 99 datapoint are females, 57 males
tapply(data$species, data$sex, function(x) length(unique(x))) # F 34; M 21 number of species
table(data$species, data$sex)
# We have thirty-eight species, with more data on females than males.
# And within the same species, there are different numbers of observations of males and females (see table(data$species, data$sex))
# Some species have only one sex (females without males).

mean_males <- mean(as.numeric(data$steepness[data$sex == "males"]), na.rm = TRUE)
mean_females <- mean(as.numeric(data$steepness[data$sex == "females"]), na.rm = TRUE)
mean_difference <- mean_males - mean_females
# In the raw data, the average steepness values are very close: males 0.81, females 0.84, difference 0.03

male_steepness <- as.numeric(data$steepness[data$sex == "males"])
female_steepness <- as.numeric(data$steepness[data$sex == "females"])
N_female_observations <- length(female_steepness)
N_male_observations   <- length(male_steepness)

plot(NA, xlim = c(0,1), ylim = c(0,7), xlab = "Steepness", ylab = "Frequency")
lines(density(male_steepness, na.rm = TRUE), col = "#FDE725FF", lwd = 8)
lines(density(female_steepness, na.rm = TRUE), col = "#443A83FF", lwd = 8)
points(rnorm(N_female_observations,mean=5,sd=0.1)~female_steepness,bg="#443A83FF",pch=21,cex=2)
points(rnorm(N_male_observations,mean=6,sd=0.1)~male_steepness,bg="#FDE725FF",pch=21,cex=2)
legend(x = "topleft", c("Males", "Females"), pch = 19, col = c("#FDE725FF", "#443A83FF"), cex = 1)
# Indeed, the difference between males and females is not striking.


# First model, straight comparison not accounting for potential dependencies among observations in the sample 
dat_list_steepness <- list(
  steepness = as.numeric(c(female_steepness, male_steepness)),  
  sex = c(rep(1, N_female_observations), rep(2, N_male_observations)))

# a very simple first model
# We assume that there is not one single mean, but two, one for each of the sexes, and determine whether these means are estimated to be different
# We need to provide priors, our expectation of what these values might be. For the means, we could expect that they are somewhere around 0.5
# For the variance, we expect this to be larger than than zero (so we use the dexp function) and larger than one
m_steepness <- ulam(
  alist(
    steepness ~ dbeta2(mean,variance),
    logit(mean) <-a[sex],
    a[sex] ~dnorm(0.5,1),   
    variance ~ dexp(10)
  ) , data=dat_list_steepness , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)
# We extract samples from the Bayesian model
post_steepness <- extract.samples(m_steepness)
# We calculated the likely means - column 1 is for the females and column 2 is for the males (because that's how we coded the data)
post_steepness$a
# [The likely variance - we assume that the variance is the same for females and for males]
post_steepness$variance
# [The model used a logit function to force the mean to be larger than zero. We now reconvert this to the actual steepness scale]
mean_females <- inv_logit(post_steepness$a[,1])
mean_males <- inv_logit(post_steepness$a[,2])
# The means we obtain here are slightly smaller than the means in the raw data (females 0.80, males 0.77) because the model takes into account that our data is not normally distributed but skewed and that values cannot be larger than 1 - but the difference between the values for the females and the males is the same (0.03)

# We calculate whether the estimated means are different or whether the distributions overlap
difference_steepness <- inv_logit(post_steepness$a[,2]) - inv_logit(post_steepness$a[,1])
results_steepness<-list(mean_females=mean_females,mean_males=mean_males,difference_steepness=difference_steepness)

precis(results_steepness)
# With this model, the estimate for the difference in the steepness of females and males crosses 0, meaning one is not consistently larger than the other : so the steepness values of the females and males are not different



#### A more complicated model that accounts for having multiple observations per species
# [We expect that steepness values for the same sex from the same species are probably similar, because they reflect the same social system]

data$sexspecies<-paste(data$species,data$sex,sep="_")

# We now combine the data into a list to be used by the model. 
# Because we have a nested model (populations are nested within species, and entries from species a grouped by whether they are from females or males),
# the entries in the list do not all come from the raw data frame, but reflect this nested structured
# That means that we designate the sex for each of the species entries rather than for each of the steepness values. We can get the number of entries from the table.
dat_list_steepness_species <- list(
  steepness = data$steepness,  # we have the steepness values in a long line, first for the females than the males
  sex = c(rep(1,tapply(data$species, data$sex, function(x) length(unique(x))) [1]),rep(2,tapply(data$species, data$sex, function(x) length(unique(x))) [2])), # we now provide the identifier that describes for each of the steepness values in the list whether it is an observation from a female (1) or a male (2) hierarchy,
  species = as.integer(as.factor(data$sexspecies))
)

# We now build our model. Again, our model goes through the steps that we used to simulate the data in reverse.
# We first assume that the steepness values come from the beta distribution with means and variances
# There is not a single mean, but a distribution that reflects that each species/sex combination can be different
# all the species-specific means for the females however should be similar, same as the species-specific means for the males, so we set is such that there are two overall means,
m_steepness_species <- ulam(
  alist(
    steepness ~ dbeta2(overallmean,amongspeciesvariance),
    logit(overallmean) <-b[species],
    b[species]~dnorm(sexspecificmean,withinspeciesvariance),
    sexspecificmean<-a[sex],
    a[sex]~dnorm(0.5,1),
    withinspeciesvariance~dexp(1),
    amongspeciesvariance~dexp(0.5)
  ) , data=dat_list_steepness_species , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)



# This is a Bayesian model, so the inference is based on sampling from the most likely space of solutions. We now work with these samples. We first extract a subset of the samples 
post_steepness_species <- extract.samples(m_steepness_species)

# We calculate the estimated mean for the females and the males. The model used a logit function to force the mean to be larger than zero. We now reconvert this to the actual steepness scale
mean_females <- inv_logit(post_steepness_species$a[,1])
mean_males <- inv_logit(post_steepness_species$a[,2])

# Next, we calculate the difference between the estimated mean for the males and the estimated mean for the females
difference_steepness <- inv_logit(post_steepness_species$a[,2]) - inv_logit(post_steepness_species$a[,1])
results_steepness<-list(mean_females=mean_females,mean_males=mean_males,difference_steepness=difference_steepness)

# We can now display the results. The inference is that, if the 5.5% - 94.5% interval for the difference does not cross zero, the steepness values of the females and males are different
precis(results_steepness)
# When accounting for biases in the sampling, that we have multiple observations from some species but not from others, the difference between females and males declines even further



##############################################################################
# Phylogenetic analyses

# We can now run the models that account for the shared phylogenetic history among species
# We first check for the phylogenetic signal, assuming that steepness values in females and in males have separate histories


mdata_phylogeny_both <- list(
  steepness_females=steepnessdata[steepnessdata$sex=="females",]$steepness,
  species_females=as.integer(as.factor((steepnessdata[steepnessdata$sex=="females",]$species))),
  N_spp_females=length(unique(steepnessdata[steepnessdata$sex=="females",]$species)),
  steepness_males=steepnessdata[steepnessdata$sex=="males",]$steepness,
  species_males=as.integer(as.factor((steepnessdata[steepnessdata$sex=="males",]$species))),
  N_spp_males=length(unique(steepnessdata[steepnessdata$sex=="males",]$species))
)

Dmat<-cophenetic(mtree)
mdata_phylogeny_both$Dmat_females<-Dmat[ unique(steepnessdata[steepnessdata$sex=="females",]$species),unique(steepnessdata[steepnessdata$sex=="females",]$species) ]/max(Dmat)
colnames(mdata_phylogeny_both$Dmat_females)<-as.integer(as.factor(colnames(mdata_phylogeny_both$Dmat_females)))
rownames(mdata_phylogeny_both$Dmat_females)<-as.integer(as.factor(rownames(mdata_phylogeny_both$Dmat_females)))
mdata_phylogeny_both$Dmat_males<-Dmat[ unique(steepnessdata[steepnessdata$sex=="males",]$species),unique(steepnessdata[steepnessdata$sex=="males",]$species) ]/max(Dmat)
colnames(mdata_phylogeny_both$Dmat_males)<-as.integer(as.factor(colnames(mdata_phylogeny_both$Dmat_males)))
rownames(mdata_phylogeny_both$Dmat_males)<-as.integer(as.factor(rownames(mdata_phylogeny_both$Dmat_males)))


m_steepness_both <- ulam(
  alist(
    steepness_males ~ dbeta2(mean_males,variance_males),
    logit(mean_males) <-a_males+b_males[species_males],
    a_males ~dnorm(0,1),
    vector[N_spp_males]:b_males~multi_normal(0,SIGMA_males),
    matrix[N_spp_males,N_spp_males]: SIGMA_males <- cov_GPL2( Dmat_males , etasq_m , rhosq_m , 0.01 ),
    etasq_m ~ half_normal(1,0.25),
    rhosq_m ~ half_normal(3,0.25),
    variance_males~dexp(10),
    steepness_females ~ dbeta2(mean_females,variance_females),
    logit(mean_females) <-a_females+b_females[species_females],
    a_females ~dnorm(0,1),
    vector[N_spp_females]:b_females~multi_normal(0,SIGMA_females),
    matrix[N_spp_females,N_spp_females]: SIGMA_females <- cov_GPL2( Dmat_females , etasq_f , rhosq_f , 0.01 ),
    etasq_f ~ half_normal(1,0.25),
    rhosq_f ~ half_normal(3,0.25),
    variance_females~dexp(10)
  ) , data=mdata_phylogeny_both , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)


samples_m_steepness_both<-extract.samples(m_steepness_both)
contrast_steepness<-as.data.frame(inv_logit(samples_m_steepness_both$a_males)-inv_logit(samples_m_steepness_both$a_females))
precis(contrast_steepness)



# We can also calculate the phylogenetic signal for the steepness values in the two sexes

specieslist_females<-as.data.frame(matrix(unique(steepnessdata[steepnessdata$sex=="females",]$species),ncol=1,nrow=length(unique(steepnessdata[steepnessdata$sex=="females",]$species))))
colnames(specieslist_females)<-"species"
rownames(specieslist_females)<-specieslist_females$species
speciesmatching_females<-name.check(phylogeny,specieslist_females)
mtree_females<-drop.tip(phylogeny,speciesmatching_females$tree_not_data)
data_female<-steepnessdata[steepnessdata$sex=="females",]
average_species_values_females<-as.data.frame(data_female %>% group_by(species) %>% summarise(meanvalue=mean(steepness)))
values_females<-average_species_values_females$meanvalue
names(values_females)<-average_species_values_females$species
phylosig(mtree_females,values_females,method="lambda",test=TRUE)
phylosig(mtree_females,values_females,method="K",test=TRUE)


specieslist_males<-as.data.frame(matrix(unique(steepnessdata[steepnessdata$sex=="males",]$species),ncol=1,nrow=length(unique(steepnessdata[steepnessdata$sex=="males",]$species))))
colnames(specieslist_males)<-"species"
rownames(specieslist_males)<-specieslist_males$species
speciesmatching_males<-name.check(phylogeny,specieslist_males)
mtree_males<-drop.tip(phylogeny,speciesmatching_males$tree_not_data)
data_male<-steepnessdata[steepnessdata$sex=="males",]
average_species_values_males<-as.data.frame(data_male %>% group_by(species) %>% summarise(meanvalue=mean(steepness)))
values_males<-average_species_values_males$meanvalue
names(values_males)<-average_species_values_males$species
phylosig(mtree_males,values_males,method="lambda",test=TRUE)
phylosig(mtree_males,values_males,method="K",test=TRUE)

# Plot the values across the phylogeny - there is generally very little variation, but there seems to be that phylogenetic pattern indicated by the phylogenetic signal

plotTree.barplot(mtree_females,values_females)
plotTree.barplot(mtree_males,values_males)
