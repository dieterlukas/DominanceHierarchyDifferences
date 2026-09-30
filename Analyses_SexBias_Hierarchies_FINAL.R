# This is code is for the analyses of the article by Spicher, Huchard, Lukas 
# The aim is to identify whether the hierarchies of males and females differ among primates

# To address this aim, we will answer the following questions:

# A) Are our measures robust?  
#   1) Are hierarchy steepness and linearity sensitive to sparseness and nb of interactions? 
#   
#   2) Are hierarchy steepness and linearity sensitive to the nb of individuals included in the hierarchy?
#
#   3) Are hierarchy steepness and linearity different in captivity than in the wild?

# B) Are there dominance types?  
#   4) Are hierarchy steepness and linearity correlated across species (w. or without accounting for phylogenetic relatedness)?
#
#   5) Are hierarchies steeper and more linear when they are based on signals rather than aggression?

# C) Are there sex-differences in the characteristics of hierarchies?
#   6) Is hierarchy steepness in males different from that in females (w. or without accounting for phylogenetic relatedness)?
#   
#   7) If steepness differs in males and females, can it be linked to the fact that female hierarchies often include more individuals?
#   
#   8) Is hierarchy linearity in males different from that in females (w. or without accounting for phylogenetic relatedness)?
#   
#   8bis) Id hierarchy differs in males and females, can it be linked to the relationship between linearity and group size ?
#
#   9) Is hierarchy linearity linked to hierarchy steepness in the same way in males and in females?
#   
#   10) Are hierarchies steeper for the sex that wins more fights - that is, is the proportion of intersexual fights that females win negatively related to the hierarchy steepness in males, and positively to the hierarchy steepness in females?
#   
#   11) Are hierarchies for females more likely to be based on signals whereas those in males more likely to be based on aggression?

# We first load the required packages
library(rethinking)
library(dplyr)
library(ape)
library(geiger)
library(phytools)
library(ggplot2)
library(gridExtra)
library(grid)
library(ggtree)
library(stringr)
library(tidyverse)
library(tidytree)
library(wesanderson)
library(cowplot)


# Load the data----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
data <- read.csv("Comparative_primate_steepness_values_clean.csv")

# Integrate the additional variable, on the intersexual dominance relationship, to the dataset (the averaged percentage of fights won by females).
dominancedata<-read.csv(url("https://github.com/dieterlukas/primate_power/raw/refs/heads/main/data/PopulationData_IntersexDominance.csv"))
intersex_dominance<-as.data.frame(dominancedata %>% group_by(species=corrected_species_id) %>% summarise(perc_won_females=mean(perc_won_females,na.rm=T)))
data<-left_join(data,intersex_dominance,by="species")
# Recover the averaged percentage of fights won by females in dominancedata for Sapajus_apella and and assign it in data for the species Cebus_apella
dominancedata %>%
  filter(corrected_species_id == "Sapajus_apella") %>%
  summarise(mean_perc = mean(perc_won_females, na.rm = TRUE))
data$perc_won_females[data$species == "Cebus_apella"] <- 33
# Some species that we have in data do not appear in dominancedata, which produces NA entries. I will delete these lines in section 10). 

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








# Analyses for section A - robustness----

# 1a) steepness sensitive to sparseness (check the sparseness for outliers first)----
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
data$studyperiod_month<-as.numeric(data$studyperiod_month)
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



# 1b) steepness sensitive to number of interactions----
dat_list_steepness_numberofineractions <- list(
  steepness = as.numeric(data$steepness),  
  numberofineractions = log(data$numberofineractions)
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


# 1c) linearity sensitive to sparseness----
dat_list_h_index_sparseness <- list(
  h_index = as.numeric(data$h_index)-0.0001,  
  sparseness = as.numeric(data$sparseness)
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


# 1d) linearity sensitive to number of interactions----
dat_list_h_index_numberofineractions <- list(
  h_index = as.numeric(data$h_index)-0.0001,  
  numberofineractions = log(data$numberofineractions)
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

# 1e) Is DCI senstive to data sparseness? ----
dat_list_dci_sparseness <- list(
  dci = as.numeric(data$dci)-0.0001,  
  sparseness = standardize(data$sparseness)
)

m_dci_sparseness <- ulam(
  alist(
    dci ~ dbeta2(mean,variance),
    logit(mean) <-a + b*sparseness,
    a ~dnorm(0.5,1),   
    b ~dnorm(0,1),
    variance ~ dexp(10)
  ) , data=dat_list_dci_sparseness , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)
# We check the results - we are interested in the effect measured in the factor b
precis(m_dci_sparseness)
plot(data$dci~data$sparseness)
#no effect of sparseness on the dci, as suggesting by Koenig et al., 2013 

# 1f) Is DCI senstive to the density of the matrice (i.e. number of interactions)? ----
dat_list_dci_density <-list(
  dci = as.numeric(data$dci)-0.0001,  
  numberofineractions = log(data$numberofineractions)
)

m_dci_density <- ulam(
  alist(
    dci ~ dbeta2(mean,variance),
    logit(mean) <-a + b*numberofineractions,
    a ~dnorm(0.5,1),   
    b ~dnorm(0,1),
    variance ~ dexp(10)
  ) , data=dat_list_dci_density , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)
# We check the results - we are interested in the effect measured in the factor b
precis(m_dci_density)
plot(data$dci~log(data$numberofineractions))
#DCI decreases as the number of recorded interactions increases
#The 89% credible interval is entirely negative. This suggests a reasonably credible effect


# 2a) steepness sensitive to the nb of individuals----
dat_list_steepness_numberofindividuals <- list(
  steepness = as.numeric(data$steepness),
  numberofindividuals = standardize(data$numberofindividuals)
)

# We need to provide priors, our expectation of what these values might be. For the means, we could expect that they are somewhere around 0.5
# For the variance, we expect this to be larger than than zero (so we use the dexp function) and larger than one
# For the relationship with numberofindividuals, we don't expect any effect a priori
m_steepness_numberofindividuals <- ulam(
  alist(
    steepness ~ dbeta2(mean,variance),
    logit(mean) <-a + b*numberofindividuals,
    a ~dnorm(0.5,1),   
    b ~dnorm(0,1),
    variance ~ dexp(10)
  ) , data=dat_list_steepness_numberofindividuals , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)
# We check the results - we are interested in the effect measured in the factor b
precis(m_steepness_numberofindividuals)
plot(data$steepness~standardize(data$numberofindividuals))
# Here, the relationship is linear, so we might want to include this as a linear predictor


# 2b) linearity sensitive to the nb of individuals----
dat_list_h_index_numberofindividuals <- list(
  h_index = as.numeric(data$h_index)-0.0001,
  numberofindividuals = standardize(data$numberofindividuals)
)

# We need to provide priors, our expectation of what these values might be. For the means, we could expect that they are somewhere around 0.5
# For the variance, we expect this to be larger than than zero (so we use the dexp function) and larger than one
# For the relationship with numberofindividuals, we don't expect any effect a priori
m_h_index_numberofindividuals <- ulam(
  alist(
    h_index ~ dbeta2(mean,variance),
    logit(mean) <-a + b*numberofindividuals,
    a ~dnorm(0.5,1),   
    b ~dnorm(0,1),
    variance ~ dexp(10)
  ) , data=dat_list_h_index_numberofindividuals , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)
# We check the results - we are interested in the effect measured in the factor b
precis(m_h_index_numberofindividuals)
plot(data$h_index~standardize(data$numberofindividuals))
# Here is no effect, so we do not need to account for it


# 2c) DCI sensitive to the nb of individuals----
dat_list_dci_numberofindividuals <- list(
  dci = as.numeric(data$dci)-0.0001,
  numberofindividuals = standardize(data$numberofindividuals)
)

# We need to provide priors, our expectation of what these values might be. For the means, we could expect that they are somewhere around 0.5
# For the variance, we expect this to be larger than than zero (so we use the dexp function) and larger than one
# For the relationship with numberofindividuals, we don't expect any effect a priori
m_dci_numberofindividuals <- ulam(
  alist(
    dci ~ dbeta2(mean,variance),
    logit(mean) <-a + b*numberofindividuals,
    a ~dnorm(0.5,1),   
    b ~dnorm(0,1),
    variance ~ dexp(10)
  ) , data=dat_list_dci_numberofindividuals , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)
# We check the results - we are interested in the effect measured in the factor b
precis(m_dci_numberofindividuals)
plot(data$dci~standardize(data$numberofindividuals))
#no effect


# 3a) steepness different in captivity than in the wild----
dat_list_steepness_captivity <- list(
  steepness = as.numeric(data$steepness),
  captive = as.numeric(as.factor(data$group)),
  numberofindividuals= standardize(data$numberofindividuals)
)

# We need to provide priors, our expectation of what these values might be. For the means, we could expect that they are somewhere around 0.5
# For the variance, we expect this to be larger than than zero (so we use the dexp function) and larger than one
# For the relationship with captivity, we don't expect any effect a priori
m_steepness_captivity <- ulam(
  alist(
    steepness ~ dbeta2(mean,variance),
    logit(mean) <-a + b[captive]+c*numberofindividuals,
    a ~dnorm(0.5,1),   
    b[captive] ~dnorm(0,1),
    c ~dnorm(0,1),
    variance ~ dexp(10)
  ) , data=dat_list_steepness_captivity , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)
# We check the results - we are interested in the effect measured in the factor b
posterior_captivity<-extract.samples(m_steepness_captivity)
contrast_captivity<-inv_logit(posterior_captivity$b[,2])-inv_logit(posterior_captivity$b[,1])
precis(contrast_captivity)


# 3b) linearity different in captivity than in the wild----
dat_list_h_index_captivity <- list(
  h_index = as.numeric(data$h_index)-0.0001,  
  captive = as.numeric(as.factor(data$group)),
  numberofindividuals= standardize(data$numberofindividuals)
)
# We need to provide priors, our expectation of what these values might be. For the means, we could expect that they are somewhere around 0.5
# For the variance, we expect this to be larger than than zero (so we use the dexp function) and larger than one
# For the relationship with captivity, we don't expect any effect a priori
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

# 3c) Is DCI different in captivity vs in wild? ----
dat_list_dci_captivity <- list(
  dci = as.numeric(data$dci)-0.0001,  
  captive = as.numeric(as.factor(data$group)))

m_dci_captivity <- ulam(
  alist(
    dci ~ dbeta2(mean,variance),
    logit(mean) <-a + b[captive],
    a ~dnorm(0.5,1),   
    b[captive] ~dnorm(0,1),
    c ~dnorm(0,1),
    variance ~ dexp(10)
  ) , data=dat_list_dci_captivity , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)
# We check the results - we are interested in the effect measured in the factor b
posterior_captivity<-extract.samples(m_dci_captivity)
contrast_captivity<-inv_logit(posterior_captivity$b[,2])-inv_logit(posterior_captivity$b[,1])
precis(contrast_captivity)






#FIGURE 2----

data2 <- data %>%
  mutate(
    steepness = as.numeric(as.character(steepness)),
    h_index = as.numeric(as.character(h_index)),
    dci = as.numeric(as.character(dci))
  )

# 1. Stats PAR CLADES (mâles+femelles)
clade_stats <- data2 %>%
  group_by(clade) %>%
  summarise(
    N_matrices = n(),
    N_species = n_distinct(species),
    
    median_steepness = round(median(steepness, na.rm = TRUE), 3),
    mean_steepness = round(mean(steepness, na.rm = TRUE), 3),
    ci_lower_steep = round(quantile(steepness, 0.10, na.rm = TRUE), 3),
    ci_upper_steep = round(quantile(steepness, 0.90, na.rm = TRUE), 3),
    
    median_linearity = round(median(h_index, na.rm = TRUE), 3),
    mean_linearity = round(mean(h_index, na.rm = TRUE), 3),
    ci_lower_lin = round(quantile(h_index, 0.10, na.rm = TRUE), 3),
    ci_upper_lin = round(quantile(h_index, 0.90, na.rm = TRUE), 3),
    
    # DCI
    median_dci = round(median(dci, na.rm=T), 3),
    mean_dci= round(mean(dci, na.rm=T), 3),
    ci_lower_steep = round(quantile(dci, 0.10, na.rm=T), 3),
    ci_upper_steep = round(quantile(dci, 0.90, na.rm=T), 3),
    
    .groups = "drop"
  )
print(clade_stats)

label_map <- c("Cercopithecoidea" = "Cercopithecoidea", 
               "Hominoidea" = "Hominoidea", 
               "Platyrrhini" = "Platyrrhini")

data2 <- data %>%
  mutate(
    steepness = as.numeric(as.character(steepness)),
    h_index = as.numeric(as.character(h_index))
  )
#plot steepness median
p_clades_steep_med <- ggplot() +
  geom_jitter(data = data2,
              aes(x = clade, y = steepness, color = clade),
              width = 0.25, size = 5, alpha = 0.6, stroke = 0) +
  # Errorbar CI médiane
  geom_errorbar(data = clade_stats,
                aes(x = clade, ymin = ci_lower_steep, ymax = ci_upper_steep),
                width = 0.25, size = 0.5, color = "black") +
  # Point médiane
  geom_point(data = clade_stats,
             aes(x = clade, y = median_steepness),
             size = 5, shape = 21, fill = "white", color = "black", stroke = 2) +
  scale_color_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb")) + 
  scale_x_discrete(labels = label_map) +
  labs(x = NULL, y = "Hierarchy Steepness") +
  ylim(0, 1) +
  theme_minimal(base_size = 16) +
  theme(axis.title = element_text(face = "bold", size = 16),
        axis.text = element_text(size = 10),
        legend.position = "none")


p_clades_lin_med <- ggplot() +
  geom_jitter(data = data2,
              aes(x = clade, y = h_index, color = clade),
              width = 0.25, size = 5, alpha = 0.6, stroke = 0) +
  # Errorbar CI médiane
  geom_errorbar(data = clade_stats,
                aes(x = clade, ymin = ci_lower_lin, ymax = ci_upper_lin),
                width = 0.25, size = 0.5, color = "black") +
  # Point médiane
  geom_point(data = clade_stats,
             aes(x = clade, y = median_linearity),
             size = 5, shape = 21, fill = "white", color = "black", stroke = 2) +
  scale_color_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb")) + 
  labs(x = "Clade", y = "Hierarchy Linearity") +
  ylim(0, 1) +
  theme_minimal(base_size = 16) +
  theme(
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 10),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

p_clades_dci_med <- ggplot() +
  geom_jitter(data = data2,
              aes(x = clade, y = dci, color = clade),
              width = 0.25, size = 5, alpha = 0.6, stroke = 0) +
  # Errorbar CI médiane
  geom_errorbar(data = clade_stats,
                aes(x = clade, ymin = ci_lower_steep, ymax = ci_upper_steep),
                width = 0.25, size = 0.5, color = "black") +
  # Point médiane
  geom_point(data = clade_stats,
             aes(x = clade, y = median_dci),
             size = 5, shape = 21, fill = "white", color = "black", stroke = 2) +
  scale_color_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb")) + 
  scale_x_discrete(labels = label_map) +
  labs(x = NULL, y = "DCI") +
  ylim(0, 1.01) +
  theme_minimal(base_size = 16) +
  theme(axis.title = element_text(face = "bold", size = 16),
        axis.text = element_text(size = 10),
        legend.position = "none")


#merge them - horizontal
figure2_horizontal <- ggdraw() +
  draw_plot(
    plot_grid(
      p_clades_steep_med + xlab(NULL),
      p_clades_lin_med + xlab(NULL),
      p_clades_dci_med + xlab(NULL),
      labels = c("A", "B", "C"),
      label_size = 16,
      label_fontface = "bold",
      ncol = 3,
      align = "hv"
    ),
    x = 0, y = 0.05, width = 1, height = 0.95
  ) +
  draw_label("Clade", x = 0.5, y = 0.03, fontface = "bold", size = 16)

figure2_horizontal

#vertical
figure2_vertical <- ggdraw() +
  draw_plot(
    plot_grid(
      p_clades_steep_med + xlab(NULL),
      p_clades_lin_med + xlab(NULL),
      p_clades_dci_med + xlab(NULL),
      labels = c("A", "B", "C"),
      ncol = 1,
      align = "v",
      label_size = 16,
      label_fontface = "bold"
    ),
    x = 0, y = 0.03, width = 1, height = 0.97
  ) +
  draw_label("Clade", x = 0.5, y = 0.01, fontface = "bold", size = 16)

figure2_vertical



#same plot but with the mean instead of median 
#plot steepness mean
p_clades_steep_mean <- ggplot() +
  geom_jitter(data = data,
              aes(x = clade, y = steepness, color = clade),
              width = 0.25, size = 5, alpha = 0.6, stroke = 0) +
  # Errorbar CI moyenne
  geom_errorbar(data = clade_stats,
                aes(x = clade, ymin = ci_lower_steep, ymax = ci_upper_steep),
                width = 0.25, size = 0.5, color = "black") +
  # Point moyen
  geom_point(data = clade_stats,
             aes(x = clade, y = mean_steepness),
             size = 5, shape = 21, fill = "white", color = "black", stroke = 2) +
  scale_color_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb")) + 
  scale_x_discrete(labels = label_map) +
  labs(x = NULL, y = "Hierarchy Steepness") +
  ylim(0, 1) +
  theme_minimal(base_size = 16) +
  theme(axis.title = element_text(face = "bold", size = 16),
        axis.text = element_text(size = 10),
        legend.position = "none")


p_clades_lin_med <- ggplot() +
  geom_jitter(data = data,
              aes(x = clade, y = h_index, color = clade),
              width = 0.25, size = 5, alpha = 0.6, stroke = 0) +
  # Errorbar CI moyenne
  geom_errorbar(data = clade_stats,
                aes(x = clade, ymin = ci_lower_lin, ymax = ci_upper_lin),
                width = 0.25, size = 0.5, color = "black") +
  # Point moyen
  geom_point(data = clade_stats,
             aes(x = clade, y = mean_linearity),
             size = 5, shape = 21, fill = "white", color = "black", stroke = 2) +
  scale_color_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb")) + 
  labs(x = "Clade", y = "Hierarchy Linearity") +
  ylim(0, 1) +
  theme_minimal(base_size = 16) +
  theme(
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 10),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

# Assemblage final
x.grob_1 <- textGrob("Clade", gp = gpar(fontsize = 16, fontface = "bold"))
grid.arrange(p_clades_steep_med + xlab(NULL), p_clades_lin_med + xlab(NULL),
             ncol = 2, bottom = x.grob_1)









# Analyses for section B - dominance types----

#   4) Are hierarchy steepness and linearity correlated across species  / accounting for phylogenetic relatedness and sparseness (which affects both)----

mdata_phylogeny_steepness_linearity <- list(
  steepness=data$steepness,
  linearity=data$h_index,
  sparseness=data$sparseness,
  species=as.integer(as.factor(data$species)),
  N_spp=length(unique(data$species))
)

mdata_phylogeny_steepness_linearity <- list(
  steepness=data$steepness,
  linearity=data$h_index-0.0001,
  sparseness=data$sparseness,
  dci=data$dci-0.0001,
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

#Are steepness and DCI correlated across species?
m_steepness_dci<- ulam(
  alist(
    steepness ~ dbeta2(mean,variance),
    logit(mean) <-a+b[species]+c*sparseness+d*dci,
    a ~dnorm(0,1),
    c~dnorm(0,1),
    d~dnorm(0,1),
    vector[N_spp]:b~multi_normal(0,SIGMA),
    matrix[N_spp,N_spp]:SIGMA <- cov_GPL2( Dmat , etasq, rhosq , 0.01 ),
    etasq ~ half_normal(1,0.25),
    rhosq ~ half_normal(3,0.25),
    variance~dexp(10)
  ) , data=mdata_phylogeny_steepness_linearity , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)

precis(m_steepness_dci)
plot(data$steepness~data$dci)
# DCI and steepness are positively correlated


#Are DCI and linearity correlated across species?
m_dci_linearity <- ulam(
  alist(
    dci ~ dbeta2(mean,variance),
    logit(mean) <-a+b[species]+c*linearity,
    a ~dnorm(0,1),
    c~dnorm(0,1),
    vector[N_spp]:b~multi_normal(0,SIGMA),
    matrix[N_spp,N_spp]:SIGMA <- cov_GPL2( Dmat , etasq, rhosq , 0.01 ),
    etasq ~ half_normal(1,0.25),
    rhosq ~ half_normal(3,0.25),
    variance~dexp(10)
  ) , data=mdata_phylogeny_steepness_linearity , chains=4 , cores=8 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)

precis(m_dci_linearity)
plot(data$dci~data$h_index)
# c crosses 0



#   5) Are hierarchies steeper and more linear when they are based on signals rather than aggression?----
data$typeofbehaviour[is.na(data$typeofbehaviour)] <- "AD"

table(data$typeofbehaviour)
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
# No real difference - difference is even smaller when accounting for the effect of sparseness on steepness


dat_list_linearity_typeofbehaviour <- list(
  linearity = as.numeric(data$h_index)-0.001,  
  typeofbehaviour = as.integer(as.factor(data$typeofbehaviour))
)
m_linearity_typeofbehaviour <- ulam(
  alist(
    linearity ~ dbeta2(mean,variance),
    logit(mean) <-a + b[typeofbehaviour],
    a ~dnorm(0.5,1),   
    b[typeofbehaviour] ~dnorm(0,1),
    variance ~ dexp(10)
  ) , data=dat_list_linearity_typeofbehaviour , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)
# We check the results - we are interested in the effect measured in the factor b
precis(m_linearity_typeofbehaviour)
post_linearity_typeofbehaviour<-extract.samples(m_linearity_typeofbehaviour)
# Comparing D versus A
contrast_D_A<-post_linearity_typeofbehaviour$b[,3]-post_linearity_typeofbehaviour$b[,1]
precis(contrast_D_A)
# No difference



dat_list_dci_typeofbehaviour <- list(
  dci = as.numeric(data$dci-0.0001),
  typeofbehaviour = as.integer(as.factor(data$typeofbehaviour))
)

m_dci_typeofbehaviour <- ulam(
  alist(
    dci ~ dbeta2(mean, variance),
    logit(mean) <- a + b[typeofbehaviour],
    a ~ dnorm(0.5, 1),
    b[typeofbehaviour] ~ dnorm(0, 1),
    variance ~ dexp(10)
  ),
  data = dat_list_dci_typeofbehaviour,
  chains = 4,
  cores = 4,
  log_lik = TRUE,
  cmdstan = TRUE,
  messages = FALSE,
  refresh = 0
)

precis(m_dci_typeofbehaviour)
post_dci_typeofbehaviour <- extract.samples(m_dci_typeofbehaviour)

contrast_D_A_dci <- post_dci_typeofbehaviour$b[,3]-post_dci_typeofbehaviour$b[,1]

precis(contrast_D_A_dci)


# Analyses for section C - sex-differences---- 
#   6) Is hierarchy steepness in males different from that in females (w. or without accounting for phylogenetic relatedness)?----

# descriptive statistics
length(unique(data$species)) # 38
table(data$sex) # 99 datapoint are females, 57 males
tapply(data$species, data$sex, function(x) length(unique(x))) # F 34; M 21 number of species
table(data$species, data$sex)
# We have thirty-eight species, with more data on females than males.
# And within the same species, there are different numbers of observations of males and females (see table(data$species, data$sex))
# Some species have only one sex (females without males).

data %>%
  group_by(sex) %>%
  summarise(
    n       = n(),
    Moyenne = round(mean(steepness, na.rm = TRUE), 3),
    SD      = round(sd(steepness, na.rm = TRUE), 3),
    Min     = round(min(steepness, na.rm = TRUE), 3),
    Max     = round(max(steepness, na.rm = TRUE), 3),
    Range   = paste0('[', Min, ' - ', Max, ']'),
    .groups = 'drop'
  )



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

#plot 2
cols_sex <- c("males"="gold3", "females"="#443A83FF")

ggplot(data, aes(x = sex, y = steepness, color = sex, fill = sex)) +
  geom_violin(alpha = 0.2, width = 0.3, color = NA) +
  stat_summary(fun = median, geom = "point", size = 4, color = "black") +
  geom_jitter(width = 0.15, size = 3, alpha = 0.6, shape = 21, stroke = 0.5) +
  scale_color_manual(values = cols_sex) +
  scale_fill_manual(values = cols_sex) +
  labs(x = "Sex", y = "Hierarchy Steepness", title = "Steepness by Sex") +
theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    legend.position = c(0, 1),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = scales::alpha("white", 0.8), color = NA)
  )

#3rd plot
height_points <- c("males" = 6, "females" = 7.5)

steep_descri <- ggplot(data, aes(x = steepness, color = sex, fill = sex)) +
  geom_density(size = 2, alpha = 0.3) +
  geom_jitter(data = data,
              aes(y = height_points[sex]),
              width = 0, size = 3, alpha = 0.6, shape = 21, fill = "white", stroke = 1.5) +
  scale_color_manual(values = cols_sex) +
  scale_fill_manual(values = cols_sex) +
  labs(x = "Hierarchy Steepness", y = "Density") +
  xlim(0,1) +
  ylim(0,8) + 
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face="bold", hjust=0.5),
    axis.title = element_text(face="bold", size=16),
    axis.text = element_text(size=14),
    legend.position="top"
  )

steep_descri

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

data$sexspecies<-paste(data$sex,data$species,sep="_")

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




### Phylogenetic analyses for 6----

# We can now run the models that account for the shared phylogenetic history among species
# We first check for the phylogenetic signal, assuming that steepness values in females and in males have separate histories
steepnessdata<-data


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

# The phylogenetic signal is stronger for the female values than for the male values. The lower value for the male values appears to occur because there is relatively little variation among species.




# [For the analysis, we now want to take into account that observations of steepness for either sex are likely to be more similar when they are from the same species. 
# We can however not simply account for species identity in this case. In each species, males and females have different social systems. Accordingly, in a given species the steepness values for females and for males can change independently. 
# Knowing, for example, that in chimpanzees steepness values are lower than the average in males does not provide any information for what the steepness values in female chimpanzees will be. 
# We therefore need a sex-specific species variable, that groups together only the observations from a single sex in a given species. We can get this by creating a new variable that combines the species name with the sex]



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


#PLOT posteriors
#transform
steepness_males_post <- inv_logit(samples_m_steepness_both$a_males)
steepness_females_post <- inv_logit(samples_m_steepness_both$a_females)

#long
plot_data <- data.frame(
  steepness = c(steepness_males_post, steepness_females_post),
  sex = rep(c("Males", "Females"), each = length(steepness_males_post))
)
 
ggplot(plot_data, aes(x = steepness, fill = sex, color = sex)) +
  geom_density(alpha = 0.3, size = 1.2) +
  scale_fill_manual(values = c("Males" = "gold3", "Females" = "#443A83FF")) +
  scale_color_manual(values = c("Males" = "gold3", "Females" = "#443A83FF")) +
  xlim(0, 1) +
  xlab("Hierarchy steepness") +
  ylab("Posterior density") +
  theme_bw(base_size = 16) +
  theme(
    legend.title = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold")
  ) 


### Checking sex differences within species----

# We can calculate the differences in steepness between the sexes within species where we have data for both females and males
steepness_withinsexcomparison<-as.data.frame(data %>% group_by(species,sex) %>% summarise(steepness=mean(steepness)))
steepness_withinsexcomparison<-steepness_withinsexcomparison %>% spread(sex,steepness)
steepness_withinsexcomparison$sexdifference<-steepness_withinsexcomparison$males-steepness_withinsexcomparison$females

# positive values mean that males have steeper hierarchies in that species, negative values females
hist(steepness_withinsexcomparison$sexdifference)
mean(steepness_withinsexcomparison$sexdifference,na.rm=T) # -0.0009515069, so not different from zero

# We can calculate all pairwise sex differences from across speciess to see whether the differences within species are different from what we would expect based on the values across species
allsexcomparisons<-NA
count<-1
for (i in 1:nrow(steepness_withinsexcomparison)){
  if(is.na(steepness_withinsexcomparison[i,]$females)){next(i)}else{
    for (j in 1:nrow(steepness_withinsexcomparison)){
      if(is.na(steepness_withinsexcomparison[j,]$males)){next(j)}else{
        allsexcomparisons[count]<-steepness_withinsexcomparison[j,]$males-steepness_withinsexcomparison[i,]$females
        count<-count+1
      }
    }
  }
}
withinsexcomparisons<-steepness_withinsexcomparison$sexdifference
withinsexcomparisons<-withinsexcomparisons[is.na(withinsexcomparisons)==F]

dat_list_sexdifference<-list(
  steepness_difference=c(allsexcomparisons,withinsexcomparisons),
  within=c(rep(0,length(allsexcomparisons)),rep(1,length(withinsexcomparisons)))
)

m_sexdifference <- ulam(
  alist(
    steepness_difference ~ dnorm(mu,sigma),
    mu <- a+b*within,
    a ~dnorm(0,1),
    b~dnorm(0,1),
    sigma~dexp(1)
  ) , data=dat_list_sexdifference , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)

precis(m_sexdifference)
# we are interested in the value of b - the value close to 0 (with 89% intervals spanning 0) means that the differences within species are indistinguishable from the sex differences between species, a negative value would have meant that within species sex differences are smaller than those between the sexes. 



# 7) If steepness differs in males and females, can it be linked to the fact that female hierarchies often include more individuals?----

# Contrary to our prediction, the more individuals there are in a group, the steeper the hierarchy
# We though also found that the steepness values of females are on average slightly larger than those of the males, so the number of individuals and the resulting number of interactions might explain the differences in steepness

male_steepness <- as.numeric(data$steepness[data$sex == "males"])
female_steepness <- as.numeric(data$steepness[data$sex == "females"])
N_female_observations <- length(female_steepness)
N_male_observations   <- length(male_steepness)
individuals = standardize(data$numberofindividuals)

dat_list_steepness_individuals <- list(
  steepness = as.numeric(c(female_steepness, male_steepness)),  
  sex = c(rep(1, N_female_observations), rep(2, N_male_observations)),
  individuals = c(individuals[data$sex=="females"],individuals[data$sex=="males"]))

# We assume that there is not one single mean, but two, one for each of the sexes, and determine whether these means are estimated to be different
# We need to provide priors, our expectation of what these values might be. For the means, we could expect that they are somewhere around 0.5
# For the variance, we expect this to be larger than than zero (so we use the dexp function) and larger than one
m_steepness_individuals <- ulam(
  alist(
    steepness ~ dbeta2(mean,variance),
    logit(mean) <-a[sex]+b*individuals,
    a[sex] ~dnorm(0.5,1), 
    b ~ dnorm(0,1),
    variance ~ dexp(10)
  ) , data=dat_list_steepness_individuals , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)
# We extract samples from the Bayesian model
post_steepness_individuals <- extract.samples(m_steepness_individuals)
# We calculated the likely means - column 1 is for the females and column 2 is for the males (because that's how we coded the data)
post_steepness_individuals$a
# [The likely variance - we assume that the variance is the same for females and for males]
post_steepness_individuals$variance
# [The model used a logit function to force the mean to be larger than zero. We now reconvert this to the actual steepness scale]
mean_females <- inv_logit(post_steepness_individuals$a[,1])
mean_males <- inv_logit(post_steepness_individuals$a[,2])
# The means we obtain here are slightly smaller than the means in the raw data (females 0.80, males 0.77) because the model takes into account that our data is not normally distributed but skewed and that values cannot be larger than 1 - but the difference between the values for the females and the males is the same (0.03)

# We calculate whether the estimated means are different or whether the distributions overlap
difference_steepness_individuals <- inv_logit(post_steepness_individuals$a[,2]) - inv_logit(post_steepness_individuals$a[,1])
results_steepness<-list(mean_females=mean_females,mean_males=mean_males,difference_steepness=difference_steepness_individuals)

precis(results_steepness)
# The difference almost disappears when taking into account that there are more females than males per group

# though it appears that the main influence is not the number of individuals per se, but the number of interactions (usually, when there are more individuals, there are more chances for interactions)
interactions<-log(data$numberofineractions)

dat_list_steepness_interactions <- list(
  steepness = as.numeric(c(female_steepness, male_steepness)),  
  sex = c(rep(1, N_female_observations), rep(2, N_male_observations)),
  interactions = c(interactions[data$sex=="females"],interactions[data$sex=="males"]))

# We assume that there is not one single mean, but two, one for each of the sexes, and determine whether these means are estimated to be different
# We need to provide priors, our expectation of what these values might be. For the means, we could expect that they are somewhere around 0.5
# For the variance, we expect this to be larger than than zero (so we use the dexp function) and larger than one
m_steepness_interactions <- ulam(
  alist(
    steepness ~ dbeta2(mean,variance),
    logit(mean) <-a[sex]+b*interactions,
    a[sex] ~dnorm(0.5,1), 
    b ~ dnorm(0,1),
    variance ~ dexp(10)
  ) , data=dat_list_steepness_interactions , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)
# We extract samples from the Bayesian model

precis(m_steepness_interactions)
#  the number of interactions has a very strong positive influence on the steepness

post_steepness_interactions <- extract.samples(m_steepness_interactions)

# [The model used a logit function to force the mean to be larger than zero. We now reconvert this to the actual steepness scale]
mean_females <- inv_logit(post_steepness_interactions$a[,1])
mean_males <- inv_logit(post_steepness_interactions$a[,2])
# The means we obtain here are slightly smaller than the means in the raw data (females 0.80, males 0.77) because the model takes into account that our data is not normally distributed but skewed and that values cannot be larger than 1 - but the difference between the values for the females and the males is the same (0.03)

# We calculate whether the estimated means are different or whether the distributions overlap
difference_steepness_individuals <- inv_logit(post_steepness_interactions$a[,2]) - inv_logit(post_steepness_interactions$a[,1])
results_steepness<-list(mean_females=mean_females,mean_males=mean_males,difference_steepness=difference_steepness_individuals)

precis(results_steepness)
# after accounting for the number of interactions, the values for males and females are very similar and the corrected estimated means are close to 0.5 because the number of interactions absorbed most of the variation in steepness


# 8) Is hierarchy linearity in males different from that in females (w. or without accounting for phylogenetic relatedness)?----
data$h_index<-data$h_index-0.001

mean_males <- mean(as.numeric(data$h_index[data$sex == "males"]), na.rm = TRUE)
mean_females <- mean(as.numeric(data$h_index[data$sex == "females"]), na.rm = TRUE)
mean_difference <- mean_males - mean_females
# In the raw data, the average h_index values are very close: males 0.81, females 0.81, difference 0.00

# Descriptives complètes
data %>%
  group_by(sex) %>%
  summarise(
    n       = n(),
    Moyenne = round(mean(h_index, na.rm = TRUE), 3),
    SD      = round(sd(h_index, na.rm = TRUE), 3),
    Min     = round(min(h_index, na.rm = TRUE), 3),
    Max     = round(max(h_index, na.rm = TRUE), 3),
    .groups = 'drop'
  )


male_h_index <- as.numeric(data$h_index[data$sex == "males"])
female_h_index <- as.numeric(data$h_index[data$sex == "females"])
N_female_observations <- length(female_h_index)
N_male_observations   <- length(male_h_index)

plot(NA, xlim = c(0,1), ylim = c(0,7), xlab = "h_index", ylab = "Frequency")
lines(density(male_h_index, na.rm = TRUE), col = "#FDE725FF", lwd = 8)
lines(density(female_h_index, na.rm = TRUE), col = "#443A83FF", lwd = 8)
points(rnorm(N_female_observations,mean=5,sd=0.1)~female_h_index,bg="#443A83FF",pch=21,cex=2)
points(rnorm(N_male_observations,mean=6,sd=0.1)~male_h_index,bg="#FDE725FF",pch=21,cex=2)
legend(x = "topleft", c("Males", "Females"), pch = 19, col = c("#FDE725FF", "#443A83FF"), cex = 1)
# Indeed, the difference between males and females is not striking.

#plot 2
cols_sex <- c("males"="gold3", "females"="#443A83FF")

ggplot(data, aes(x = sex, y = h_index, color = sex, fill = sex)) +
  geom_violin(alpha = 0.2, width = 0.3, color = NA) +
  stat_summary(fun = median, geom = "point", size = 4, color = "black") +
  geom_jitter(width = 0.15, size = 3, alpha = 0.6, shape = 21, stroke = 0.5) +
  scale_color_manual(values = cols_sex) +
  scale_fill_manual(values = cols_sex) +
  labs(x = "Sex", y = "Hierarchy Linearity") +
  ylim(0, 1) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none"
  )


              
#3rd plot
height_points <- c("males" = 4.5, "females" = 5.5)

linearity_descri <- ggplot(data, aes(x = h_index, color = sex, fill = sex)) +
  geom_density(size = 2, alpha = 0.3) +
  geom_jitter(data = data,
              aes(y = height_points[sex]),
              width = 0, size = 3, alpha = 0.6, shape = 21, fill = "white", stroke = 1.5) +
  scale_color_manual(values = cols_sex) +
  scale_fill_manual(values = cols_sex) +
  labs(x = "Hierarchy Linearity", y = "Density") +
  xlim(0,1) +
  ylim(0,6) + 
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face="bold", hjust=0.5),
    axis.title = element_text(face="bold", size=16),
    axis.text = element_text(size=14),
    legend.position="none"
  )  
linearity_descri



#FIGURE 5----
#second trial
height_points <- c("males" = 6, "females" = 7.5)

steep_descri <- ggplot(data, aes(x = steepness, color = sex, fill = sex)) +
  geom_density(size = 2, alpha = 0.3) +
  geom_jitter(aes(y = height_points[sex]),
              width = 0, size = 3, alpha = 0.6, shape = 21,
              fill = "white", stroke = 1.5) +
  scale_color_manual(values = cols_sex) +
  scale_fill_manual(values = cols_sex) +
  labs(x = "Hierarchy Steepness", y = "Density") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 8)) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 14),
    legend.position = "top"
  )


linearity_descri <- ggplot(data, aes(x = h_index, color = sex, fill = sex)) +
  geom_density(size = 2, alpha = 0.3) +
  geom_jitter(data = data,
              aes(y = height_points[sex]),
              width = 0, size = 3, alpha = 0.6, shape = 21, fill = "white", stroke = 1.5) +
  scale_color_manual(values = cols_sex) +
  scale_fill_manual(values = cols_sex) +
  labs(x = "Hierarchy Linearity", y = "Density") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 8)) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 14),
    legend.position = "top"
  )

dci_descri <- ggplot(data, aes(x = dci, color = sex, fill = sex)) +
  geom_density(size = 2, alpha = 0.3) +
  geom_jitter(data = data,
              aes(y = height_points[sex]),
              width = 0, size = 3, alpha = 0.6, shape = 21, fill = "white", stroke = 1.5) +
  scale_color_manual(values = cols_sex) +
  scale_fill_manual(values = cols_sex) +
  labs(x = "DCI", y = "Density") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 8)) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 14),
    legend.position = "top"
  )


p1 <- steep_descri + theme(legend.position = "none")
p2 <- linearity_descri + theme(legend.position = "none")
p3 <- dci_descri + theme(legend.position = "none")

# 2) récupérer une légende commune
legend_top <- get_legend(
  steep_descri +
    guides(
      color = guide_legend(
        nrow = 1,
        override.aes = list(
          shape = 21,
          size = 4,
          fill = "white",
          alpha = 1,
          stroke = 1.5,
          linetype = 0
        )
      ),
      fill = "none"
    ) +
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.title = element_blank(),
      legend.box = "horizontal"
    )
)

# 3) aligner les 3 plots en ligne
plots_row <- plot_grid(
  p1, p2, p3,
  nrow = 1,
  align = "hv",
  labels = c("A", "B", "C"),
  label_size = 16,
  label_fontface = "bold"
)

# 4) ajouter la légende en haut
final_plot <- plot_grid(
  legend_top,
  plots_row,
  ncol = 1,
  rel_heights = c(0.12, 1)
)

final_plot

#vertical
plots_col <- plot_grid(
  p1, p2, p3,
  ncol = 1,
  align = "v",
  labels = c("A", "B", "C"),
  label_size = 16,
  label_fontface = "bold"
)

final_plot <- plot_grid(
  legend_top,
  plots_col,
  ncol = 1,
  rel_heights = c(0.12, 1)
)

final_plot



# First model, straight comparison not accounting for potential dependencies among observations in the sample 
dat_list_h_index <- list(
  h_index = as.numeric(c(female_h_index, male_h_index)),  
  sex = c(rep(1, N_female_observations), rep(2, N_male_observations)))

# a very simple first model
# We assume that there is not one single mean, but two, one for each of the sexes, and determine whether these means are estimated to be different
# We need to provide priors, our expectation of what these values might be. For the means, we could expect that they are somewhere around 0.5
# For the variance, we expect this to be larger than than zero (so we use the dexp function) and larger than one
m_h_index <- ulam(
  alist(
    h_index ~ dbeta2(mean,variance),
    logit(mean) <-a[sex],
    a[sex] ~dnorm(0.5,1),   
    variance ~ dexp(10)
  ) , data=dat_list_h_index , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)
# We extract samples from the Bayesian model
post_h_index <- extract.samples(m_h_index)
# We calculated the likely means - column 1 is for the females and column 2 is for the males (because that's how we coded the data)
post_h_index$a
# [The likely variance - we assume that the variance is the same for females and for males]
post_h_index$variance
# [The model used a logit function to force the mean to be larger than zero. We now reconvert this to the actual h_index scale]
mean_females <- inv_logit(post_h_index$a[,1])
mean_males <- inv_logit(post_h_index$a[,2])
# The means we obtain here are slightly smaller than the means in the raw data (females 0.80, males 0.77) because the model takes into account that our data is not normally distributed but skewed and that values cannot be larger than 1 - but the difference between the values for the females and the males is the same (0.03)

# We calculate whether the estimated means are different or whether the distributions overlap
difference_h_index <- inv_logit(post_h_index$a[,2]) - inv_logit(post_h_index$a[,1])
results_h_index<-list(mean_females=mean_females,mean_males=mean_males,difference_h_index=difference_h_index)

precis(results_h_index)
# With this model, the estimate for the difference in the h_index of females and males crosses 0, meaning one is not consistently larger than the other : so the h_index values of the females and males are not different



#### A more complicated model that accounts for having multiple observations per species
# [We expect that h_index values for the same sex from the same species are probably similar, because they reflect the same social system]

data$sexspecies<-paste(data$sex,data$species,sep="_")

# We now combine the data into a list to be used by the model. 
# Because we have a nested model (populations are nested within species, and entries from species a grouped by whether they are from females or males),
# the entries in the list do not all come from the raw data frame, but reflect this nested structured
# That means that we designate the sex for each of the species entries rather than for each of the h_index values. We can get the number of entries from the table.
dat_list_h_index_species <- list(
  h_index = data$h_index,  # we have the h_index values in a long line, first for the females than the males
  sex = c(rep(1,tapply(data$species, data$sex, function(x) length(unique(x))) [1]),rep(2,tapply(data$species, data$sex, function(x) length(unique(x))) [2])), # we now provide the identifier that describes for each of the h_index values in the list whether it is an observation from a female (1) or a male (2) hierarchy,
  species = as.integer(as.factor(data$sexspecies))
)

# We now build our model. Again, our model goes through the steps that we used to simulate the data in reverse.
# We first assume that the h_index values come from the beta distribution with means and variances
# There is not a single mean, but a distribution that reflects that each species/sex combination can be different
# all the species-specific means for the females however should be similar, same as the species-specific means for the males, so we set is such that there are two overall means,
m_h_index_species <- ulam(
  alist(
    h_index ~ dbeta2(overallmean,amongspeciesvariance),
    logit(overallmean) <-b[species],
    b[species]~dnorm(sexspecificmean,withinspeciesvariance),
    sexspecificmean<-a[sex],
    a[sex]~dnorm(0.5,1),
    withinspeciesvariance~dexp(1),
    amongspeciesvariance~dexp(0.5)
  ) , data=dat_list_h_index_species , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)



# This is a Bayesian model, so the inference is based on sampling from the most likely space of solutions. We now work with these samples. We first extract a subset of the samples 
post_h_index_species <- extract.samples(m_h_index_species)

# We calculate the estimated mean for the females and the males. The model used a logit function to force the mean to be larger than zero. We now reconvert this to the actual h_index scale
mean_females <- inv_logit(post_h_index_species$a[,1])
mean_males <- inv_logit(post_h_index_species$a[,2])

# Next, we calculate the difference between the estimated mean for the males and the estimated mean for the females
difference_h_index <- inv_logit(post_h_index_species$a[,2]) - inv_logit(post_h_index_species$a[,1])
results_h_index<-list(mean_females=mean_females,mean_males=mean_males,difference_h_index=difference_h_index)

# We can now display the results. The inference is that, if the 5.5% - 94.5% interval for the difference does not cross zero, the h_index values of the females and males are different
precis(results_h_index)
# When accounting for biases in the sampling, that we have multiple observations from some species but not from others, there appears to be a difference in the linearity: females have more linear hierarchies than males (note, this is still without accounting for sparseness)




### Phylogenetic analyses fir 8)----

# We can now run the models that account for the shared phylogenetic history among species
# We first check for the phylogenetic signal, assuming that h_index values in females and in males have separate histories
h_indexdata<-data


specieslist_females<-as.data.frame(matrix(unique(h_indexdata[h_indexdata$sex=="females",]$species),ncol=1,nrow=length(unique(h_indexdata[h_indexdata$sex=="females",]$species))))
colnames(specieslist_females)<-"species"
rownames(specieslist_females)<-specieslist_females$species
speciesmatching_females<-name.check(phylogeny,specieslist_females)
mtree_females<-drop.tip(phylogeny,speciesmatching_females$tree_not_data)
data_female<-h_indexdata[h_indexdata$sex=="females",]
average_species_values_females<-as.data.frame(data_female %>% group_by(species) %>% summarise(meanvalue=mean(h_index)))
values_females<-average_species_values_females$meanvalue
names(values_females)<-average_species_values_females$species
phylosig(mtree_females,values_females,method="lambda",test=TRUE)
phylosig(mtree_females,values_females,method="K",test=TRUE)


specieslist_males<-as.data.frame(matrix(unique(h_indexdata[h_indexdata$sex=="males",]$species),ncol=1,nrow=length(unique(h_indexdata[h_indexdata$sex=="males",]$species))))
colnames(specieslist_males)<-"species"
rownames(specieslist_males)<-specieslist_males$species
speciesmatching_males<-name.check(phylogeny,specieslist_males)
mtree_males<-drop.tip(phylogeny,speciesmatching_males$tree_not_data)
data_male<-h_indexdata[h_indexdata$sex=="males",]
average_species_values_males<-as.data.frame(data_male %>% group_by(species) %>% summarise(meanvalue=mean(h_index)))
values_males<-average_species_values_males$meanvalue
names(values_males)<-average_species_values_males$species
phylosig(mtree_males,values_males,method="lambda",test=TRUE)
phylosig(mtree_males,values_males,method="K",test=TRUE)

# Plot the values across the phylogeny - there is generally very little variation, but there seems to be that phylogenetic pattern indicated by the phylogenetic signal

plotTree.barplot(mtree_females,values_females)
plotTree.barplot(mtree_males,values_males)

# For linearity, there is no real phylogenetic signal, what little there is is stronger in males than in females. In this case, the signal appears absent not because there is no variatio in linearity, but because closely related species can differ as much as more distantly related species.



# [For the analysis, we now want to take into account that observations of h_index for either sex are likely to be more similar when they are from the same species. 
# We can however not simply account for species identity in this case. In each species, males and females have different social systems. Accordingly, in a given species the h_index values for females and for males can change independently. 
# Knowing, for example, that in chimpanzees h_index values are lower than the average in males does not provide any information for what the h_index values in female chimpanzees will be. 
# We therefore need a sex-specific species variable, that groups together only the observations from a single sex in a given species. We can get this by creating a new variable that combines the species name with the sex]



mdata_phylogeny_both <- list(
  h_index_females=h_indexdata[h_indexdata$sex=="females",]$h_index,
  species_females=as.integer(as.factor((h_indexdata[h_indexdata$sex=="females",]$species))),
  N_spp_females=length(unique(h_indexdata[h_indexdata$sex=="females",]$species)),
  h_index_males=h_indexdata[h_indexdata$sex=="males",]$h_index,
  species_males=as.integer(as.factor((h_indexdata[h_indexdata$sex=="males",]$species))),
  N_spp_males=length(unique(h_indexdata[h_indexdata$sex=="males",]$species))
)

Dmat<-cophenetic(mtree)
mdata_phylogeny_both$Dmat_females<-Dmat[ unique(h_indexdata[h_indexdata$sex=="females",]$species),unique(h_indexdata[h_indexdata$sex=="females",]$species) ]/max(Dmat)
colnames(mdata_phylogeny_both$Dmat_females)<-as.integer(as.factor(colnames(mdata_phylogeny_both$Dmat_females)))
rownames(mdata_phylogeny_both$Dmat_females)<-as.integer(as.factor(rownames(mdata_phylogeny_both$Dmat_females)))
mdata_phylogeny_both$Dmat_males<-Dmat[ unique(h_indexdata[h_indexdata$sex=="males",]$species),unique(h_indexdata[h_indexdata$sex=="males",]$species) ]/max(Dmat)
colnames(mdata_phylogeny_both$Dmat_males)<-as.integer(as.factor(colnames(mdata_phylogeny_both$Dmat_males)))
rownames(mdata_phylogeny_both$Dmat_males)<-as.integer(as.factor(rownames(mdata_phylogeny_both$Dmat_males)))


m_h_index_both <- ulam(
  alist(
    h_index_males ~ dbeta2(mean_males,variance_males),
    logit(mean_males) <-a_males+b_males[species_males],
    a_males ~dnorm(0,1),
    vector[N_spp_males]:b_males~multi_normal(0,SIGMA_males),
    matrix[N_spp_males,N_spp_males]: SIGMA_males <- cov_GPL2( Dmat_males , etasq_m , rhosq_m , 0.01 ),
    etasq_m ~ half_normal(1,0.25),
    rhosq_m ~ half_normal(3,0.25),
    variance_males~dexp(10),
    h_index_females ~ dbeta2(mean_females,variance_females),
    logit(mean_females) <-a_females+b_females[species_females],
    a_females ~dnorm(0,1),
    vector[N_spp_females]:b_females~multi_normal(0,SIGMA_females),
    matrix[N_spp_females,N_spp_females]: SIGMA_females <- cov_GPL2( Dmat_females , etasq_f , rhosq_f , 0.01 ),
    etasq_f ~ half_normal(1,0.25),
    rhosq_f ~ half_normal(3,0.25),
    variance_females~dexp(10)
  ) , data=mdata_phylogeny_both , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)


samples_m_h_index_both<-extract.samples(m_h_index_both)
contrast_h_index<-as.data.frame(inv_logit(samples_m_h_index_both$a_males)-inv_logit(samples_m_h_index_both$a_females))
precis(contrast_h_index)




#### A more complicated model that accounts for having multiple observations per species
# [We expect that h_index values for the same sex from the same species are probably similar, because they reflect the same social system] 
# and accounting that sparseness affects the estimate

data$sexspecies<-paste(data$sex,data$species,sep="_")

# We now combine the data into a list to be used by the model. 
# Because we have a nested model (populations are nested within species, and entries from species a grouped by whether they are from females or males),
# the entries in the list do not all come from the raw data frame, but reflect this nested structured
# That means that we designate the sex for each of the species entries rather than for each of the h_index values. We can get the number of entries from the table.
dat_list_h_index_species <- list(
  h_index = data$h_index,  # we have the h_index values in a long line, first for the females than the males
  sex = c(rep(1,tapply(data$species, data$sex, function(x) length(unique(x))) [1]),rep(2,tapply(data$species, data$sex, function(x) length(unique(x))) [2])), # we now provide the identifier that describes for each of the h_index values in the list whether it is an observation from a female (1) or a male (2) hierarchy,
  species = as.integer(as.factor(data$sexspecies)),
  sparseness = data$sparseness
)

# We now build our model. Again, our model goes through the steps that we used to simulate the data in reverse.
# We first assume that the h_index values come from the beta distribution with means and variances
# There is not a single mean, but a distribution that reflects that each species/sex combination can be different
# all the species-specific means for the females however should be similar, same as the species-specific means for the males, so we set is such that there are two overall means,
m_h_index_species <- ulam(
  alist(
    h_index ~ dbeta2(overallmean,amongspeciesvariance),
    logit(overallmean) <-b[species]+c*sparseness,
    b[species]~dnorm(sexspecificmean,withinspeciesvariance),
    sexspecificmean<-a[sex],
    a[sex]~dnorm(0.5,1),
    c~dnorm(0,1),
    withinspeciesvariance~dexp(1),
    amongspeciesvariance~dexp(0.5)
  ) , data=dat_list_h_index_species , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)


# This is a Bayesian model, so the inference is based on sampling from the most likely space of solutions. We now work with these samples. We first extract a subset of the samples 
post_h_index_species <- extract.samples(m_h_index_species)

# We calculate the estimated mean for the females and the males. The model used a logit function to force the mean to be larger than zero. We now reconvert this to the actual h_index scale
mean_females <- inv_logit(post_h_index_species$a[,1])
mean_males <- inv_logit(post_h_index_species$a[,2])

# Next, we calculate the difference between the estimated mean for the males and the estimated mean for the females
difference_h_index <- inv_logit(post_h_index_species$a[,2]) - inv_logit(post_h_index_species$a[,1])
results_h_index<-list(mean_females=mean_females,mean_males=mean_males,difference_h_index=difference_h_index)

# We can now display the results. The inference is that, if the 5.5% - 94.5% interval for the difference does not cross zero, the h_index values of the females and males are different
precis(results_h_index)

                                                                                     
#PLOT posteriors
#transform
linearity_males_post <- inv_logit(samples_m_h_index_both$a_males)
linearity_females_post <- inv_logit(samples_m_h_index_both$a_females)

#long
plot_data <- data.frame(
  linearity = c(linearity_males_post, linearity_females_post),
  sex = rep(c("Males", "Females"), each = length(linearity_males_post))
)

ggplot(plot_data, aes(x = linearity, fill = sex, color = sex)) +
  geom_density(alpha = 0.3, size = 1.2) +
  scale_fill_manual(values = c("Males" = "gold3", "Females" = "#443A83FF")) +
  scale_color_manual(values = c("Males" = "gold3", "Females" = "#443A83FF")) +
  xlim(0, 1) +
  xlab("Hierarchy linearity") +
  ylab("Posterior density") +
  theme_bw(base_size = 16) +
  theme(
    legend.title = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold")
  ) 


# 8bis) If linearity differs in males and females, can it be linked to the relationship between linearity and group size?----
data$h_index<-data$h_index-0.001

male_linearity <- as.numeric(data$h_index[data$sex == "males"])
female_linearity <- as.numeric(data$h_index[data$sex == "females"])

N_female_observations <- length(female_linearity)
N_male_observations   <- length(male_linearity)

individuals <- standardize(data$numberofindividuals)

dat_list_linearity_individuals <- list(
  linearity = as.numeric(c(female_linearity, male_linearity)),
  sex = c(rep(1, N_female_observations), rep(2, N_male_observations)),
  individuals = c(individuals[data$sex == "females"],
                  individuals[data$sex == "males"])
)

m_linearity_individuals <- ulam(
  alist(
    linearity ~ dbeta2(mean, variance),
    logit(mean) <- a[sex] + b * individuals,
    a[sex] ~ dnorm(0.5, 1),
    b ~ dnorm(0, 1),
    variance ~ dexp(10)
  ),
  data = dat_list_linearity_individuals,
  chains = 4, cores = 4, log_lik = TRUE,
  cmdstan = TRUE, messages = FALSE, refresh = 0
)

post_linearity_individuals <- extract.samples(m_linearity_individuals)

mean_females <- inv_logit(post_linearity_individuals$a[,1])
mean_males   <- inv_logit(post_linearity_individuals$a[,2])

difference_linearity <- mean_males - mean_females
difference_linearity_individuals <- inv_logit(post_linearity_individuals$a[,2]) - inv_logit(post_linearity_individuals$a[,1])

results_linearity<-list(mean_females=mean_females,mean_males=mean_males,difference_linearity=difference_linearity_individuals)

precis(results_linearity)



# 9) Is DCI in males different from that in females (w. or without accounting for phylogenetic relatedness)?----
data %>%
  group_by(sex) %>%
  summarise(
    n       = n(),
    Moyenne = round(mean(dci, na.rm = TRUE), 3),
    SD      = round(sd(dci, na.rm = TRUE), 3),
    Min     = round(min(dci, na.rm = TRUE), 3),
    Max     = round(max(dci, na.rm = TRUE), 3),
    Range   = paste0('[', Min, ' - ', Max, ']'),
    .groups = 'drop'
  )

mean_males <- mean(as.numeric(data$dci[data$sex == "males"]), na.rm = TRUE)
mean_females <- mean(as.numeric(data$dci[data$sex == "females"]), na.rm = TRUE)
mean_difference <- mean_males - mean_females

male_dci <- as.numeric(data$dci[data$sex == "males"])
female_dci <- as.numeric(data$dci[data$sex == "females"])
N_female_observations <- length(female_dci)
N_male_observations   <- length(male_dci)

plot(NA, xlim = c(0,1), ylim = c(0,7), xlab = "DCI", ylab = "Frequency")
lines(density(male_dci, na.rm = TRUE), col = "#FDE725FF", lwd = 8)
lines(density(female_dci, na.rm = TRUE), col = "#443A83FF", lwd = 8)
points(rnorm(N_female_observations,mean=5,sd=0.1)~female_dci,bg="#443A83FF",pch=21,cex=2)
points(rnorm(N_male_observations,mean=6,sd=0.1)~male_dci,bg="#FDE725FF",pch=21,cex=2)
legend(x = "topleft", c("Males", "Females"), pch = 19, col = c("#FDE725FF", "#443A83FF"), cex = 1)

data$sexspecies<-paste(data$sex,data$species,sep="_")

dat_list_dci_species <- list(
  dci = as.numeric(data$dci)-0.0001, 
  sex = c(rep(1,tapply(data$species, data$sex, function(x) length(unique(x))) [1]),rep(2,tapply(data$species, data$sex, function(x) length(unique(x))) [2])), # we now provide the identifier that describes for each of the steepness values in the list whether it is an observation from a female (1) or a male (2) hierarchy,
  species = as.integer(as.factor(data$sexspecies))
)

# We now build our model. Again, our model goes through the steps that we used to simulate the data in reverse.
# We first assume that the steepness values come from the beta distribution with means and variances
# There is not a single mean, but a distribution that reflects that each species/sex combination can be different
# all the species-specific means for the females however should be similar, same as the species-specific means for the males, so we set is such that there are two overall means,
m_dci_species <- ulam(
  alist(
    dci ~ dbeta2(overallmean,amongspeciesvariance),
    logit(overallmean) <-b[species],
    b[species]~dnorm(sexspecificmean,withinspeciesvariance),
    sexspecificmean<-a[sex],
    a[sex]~dnorm(0.5,1),
    withinspeciesvariance~dexp(1),
    amongspeciesvariance~dexp(0.5)
  ) , data=dat_list_dci_species , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)

post_dci_species <- extract.samples(m_dci_species)

# We calculate the estimated mean for the females and the males. The model used a logit function to force the mean to be larger than zero. We now reconvert this to the actual steepness scale
mean_females <- inv_logit(post_dci_species$a[,1])
mean_males <- inv_logit(post_dci_species$a[,2])

# Next, we calculate the difference between the estimated mean for the males and the estimated mean for the females
difference_dci <- inv_logit(post_dci_species$a[,2]) - inv_logit(post_dci_species$a[,1])
results_dci<-list(mean_females=mean_females,mean_males=mean_males,difference_dci=difference_dci)

# We can now display the results. The inference is that, if the 5.5% - 94.5% interval for the difference does not cross zero, the steepness values of the females and males are different
precis(results_dci)
# When accounting for biases in the sampling, that we hav



### Phylogenetic analyses for 9)----

# We can now run the models that account for the shared phylogenetic history among species
# We first check for the phylogenetic signal, assuming that steepness values in females and in males have separate histories
dcidata<-data


specieslist_females<-as.data.frame(matrix(unique(dcidata[dcidata$sex=="females",]$species),ncol=1,nrow=length(unique(dcidata[dcidata$sex=="females",]$species))))
colnames(specieslist_females)<-"species"
rownames(specieslist_females)<-specieslist_females$species
speciesmatching_females<-name.check(phylogeny,specieslist_females)
mtree_females<-drop.tip(phylogeny,speciesmatching_females$tree_not_data)
data_female<-dcidata[dcidata$sex=="females",]
average_species_values_females<-as.data.frame(data_female %>% group_by(species) %>% summarise(meanvalue=mean(dci)))
values_females<-average_species_values_females$meanvalue
names(values_females)<-average_species_values_females$species
phylosig(mtree_females,values_females,method="lambda",test=TRUE)
phylosig(mtree_females,values_females,method="K",test=TRUE)
#Female DCI showed a moderate phylogenetic tendency, with a relatively high Pagel’s lambda but a low Blomberg’s K. 
#However, the lambda test was marginal (P = 0.051) and K was not significant, suggesting that the phylogenetic signal is present but not especially strong under a Brownian expectation



specieslist_males<-as.data.frame(matrix(unique(dcidata[dcidata$sex=="males",]$species),ncol=1,nrow=length(unique(dcidata[dcidata$sex=="males",]$species))))
colnames(specieslist_males)<-"species"
rownames(specieslist_males)<-specieslist_males$species
speciesmatching_males<-name.check(phylogeny,specieslist_males)
mtree_males<-drop.tip(phylogeny,speciesmatching_males$tree_not_data)
data_male<-dcidata[dcidata$sex=="males",]
average_species_values_males<-as.data.frame(data_male %>% group_by(species) %>% summarise(meanvalue=mean(dci)))
values_males<-average_species_values_males$meanvalue
names(values_males)<-average_species_values_males$species
phylosig(mtree_males,values_males,method="lambda",test=TRUE)
phylosig(mtree_males,values_males,method="K",test=TRUE)
#Male DCI showed a weak and statistically uncertain phylogenetic signal. 
#Pagel’s lambda was close to 1 but non-significant, while Blomberg’s K was moderate and marginally non-significant, suggesting only limited evidence that closely related species resemble each other in male DCI

# Plot the values across the phylogeny - there is generally very little variation, but there seems to be that phylogenetic pattern indicated by the phylogenetic signal

plotTree.barplot(mtree_females,values_females)
plotTree.barplot(mtree_males,values_males)

# The phylogenetic signal is stronger for the female values than for the male values. The lower value for the male values appears to occur because there is relatively little variation among species.

dci_females <- as.numeric(dcidata[dcidata$sex == "females", ]$dci) - 0.0001
dci_males   <- as.numeric(dcidata[dcidata$sex == "males", ]$dci) - 0.0001


mdata_phylogeny_both <- list(
  dci_females = dci_females,
  species_females = as.integer(as.factor(dcidata[dcidata$sex == "females", ]$species)),
  N_spp_females = length(unique(dcidata[dcidata$sex == "females", ]$species)),
  dci_males = dci_males,
  species_males = as.integer(as.factor(dcidata[dcidata$sex == "males", ]$species)),
  N_spp_males = length(unique(dcidata[dcidata$sex == "males", ]$species))
)


Dmat<-cophenetic(mtree)
mdata_phylogeny_both$Dmat_females<-Dmat[ unique(dcidata[dcidata$sex=="females",]$species),unique(dcidata[dcidata$sex=="females",]$species) ]/max(Dmat)
colnames(mdata_phylogeny_both$Dmat_females)<-as.integer(as.factor(colnames(mdata_phylogeny_both$Dmat_females)))
rownames(mdata_phylogeny_both$Dmat_females)<-as.integer(as.factor(rownames(mdata_phylogeny_both$Dmat_females)))
mdata_phylogeny_both$Dmat_males<-Dmat[ unique(dcidata[dcidata$sex=="males",]$species),unique(dcidata[dcidata$sex=="males",]$species) ]/max(Dmat)
colnames(mdata_phylogeny_both$Dmat_males)<-as.integer(as.factor(colnames(mdata_phylogeny_both$Dmat_males)))
rownames(mdata_phylogeny_both$Dmat_males)<-as.integer(as.factor(rownames(mdata_phylogeny_both$Dmat_males)))


m_dci_both <- ulam(
  alist(
    dci_males ~ dbeta2(mean_males,variance_males),
    logit(mean_males) <-a_males+b_males[species_males],
    a_males ~dnorm(0,1),
    vector[N_spp_males]:b_males~multi_normal(0,SIGMA_males),
    matrix[N_spp_males,N_spp_males]: SIGMA_males <- cov_GPL2( Dmat_males , etasq_m , rhosq_m , 0.01 ),
    etasq_m ~ half_normal(1,0.25),
    rhosq_m ~ half_normal(3,0.25),
    variance_males~dexp(10),
    dci_females ~ dbeta2(mean_females,variance_females),
    logit(mean_females) <-a_females+b_females[species_females],
    a_females ~dnorm(0,1),
    vector[N_spp_females]:b_females~multi_normal(0,SIGMA_females),
    matrix[N_spp_females,N_spp_females]: SIGMA_females <- cov_GPL2( Dmat_females , etasq_f , rhosq_f , 0.01 ),
    etasq_f ~ half_normal(1,0.25),
    rhosq_f ~ half_normal(3,0.25),
    variance_females~dexp(10)
  ) , data=mdata_phylogeny_both , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)


samples_m_dci_both<-extract.samples(m_dci_both)
contrast_dci<-as.data.frame(inv_logit(samples_m_dci_both$a_males)-inv_logit(samples_m_dci_both$a_females))
precis(contrast_dci)


#PLOT posteriors
#transform
dci_males_post <- inv_logit(samples_m_dci_both$a_males)
dci_females_post <- inv_logit(samples_m_dci_both$a_females)

#long
plot_data <- data.frame(
  dci = c(dci_males_post, dci_females_post),
  sex = rep(c("Males", "Females"), each = length(dci_males_post))
)

ggplot(plot_data, aes(x = dci, fill = sex, color = sex)) +
  geom_density(alpha = 0.3, size = 1.2) +
  scale_fill_manual(values = c("Males" = "gold3", "Females" = "#443A83FF")) +
scale_color_manual(values = c("Males" = "gold3", "Females" = "#443A83FF")) +
  xlim(0, 1) +
  xlab("DCI") +
  ylab("Posterior") +
  theme_bw(base_size = 16) +
  theme(
    legend.title = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold")
  ) 


#   10) Is hierarchy linearity linked to hierarchy steepness in the same way in males and in females?----

mdata_sex_steepness_linearity <- list(
  femalesteepness=as.numeric(data[data$sex=="females",]$steepness),
  malesteepness=as.numeric(data[data$sex=="males",]$steepness),
  femalelinearity=standardize(data[data$sex=="females",]$h_index),
  malelinearity=standardize(data[data$sex=="males",]$h_index),
  femalesparseness=standardize(data[data$sex=="females",]$sparseness),
  malesparseness=standardize(data[data$sex=="males",]$sparseness)
)


m_sex_steepness_linearity <- ulam(
  alist(
    femalesteepness ~ dbeta2(femalemean,femalevariance),
    logit(femalemean) <-af+b*femalesparseness+cf*femalelinearity,
    malesteepness ~ dbeta2(malemean,malevariance),
    logit(malemean) <-am+b*malesparseness+(cf+maleoffset)*malelinearity,
    c(af,am)~dnorm(0,1),
    b~dnorm(0,1),
    cf~dnorm(0,1),
    maleoffset~dnorm(0,1),
    femalevariance~dexp(10),
    malevariance~dexp(10)
  ) , data=mdata_sex_steepness_linearity , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)

precis(m_sex_steepness_linearity)
# the estimate for the male offset is not different from zero, suggesting that the relationship between linearity and steepness is identical in males and females

plot(data[data$sex=="females",]$steepness~data[data$sex=="females",]$h_index,col="purple")
points(data[data$sex=="males",]$steepness~data[data$sex=="males",]$h_index,col="darkgreen",pch=16)
# Linearity and steepness are positively correlated
                                                                                       
ggplot(data, aes(x = h_index, y = steepness, color = sex)) +
  geom_point(size = 3, alpha = 0.6) +
  scale_color_manual(values=c("females"="#443A83FF", "males"="gold3")) +
  xlab("Linearity (h_index)") +
  ylab("Steepness") +
  theme_minimal(base_size = 16) +
  theme(legend.position="top")

                           
#   11) Are hierarchies steeper for the sex that wins more fights - that is, is the proportion of intersexual fights that females win negatively related to the hierarchy steepness in males, and positively to the hierarchy steepness in females?----

data_dominance<-data[is.na(data$perc_won_females)==F,]
data_dominance$perc_won_females<-log(data_dominance$perc_won_females+0.01)

mdata_sex_steepness_dominance <- list(
  femalesteepness=as.numeric(data_dominance[data_dominance$sex=="females",]$steepness),
  malesteepness=as.numeric(data_dominance[data_dominance$sex=="males",]$steepness),
  dominance_femalevalues=standardize(data_dominance[data_dominance$sex=="females",]$perc_won_females),
  dominance_malevalues=standardize(data_dominance[data_dominance$sex=="males",]$perc_won_females),
  femalesparseness=standardize(data_dominance[data_dominance$sex=="females",]$sparseness),
  malesparseness=standardize(data_dominance[data_dominance$sex=="males",]$sparseness)
)

m_sex_steepness_dominance <- ulam(
  alist(
    femalesteepness ~ dbeta2(femalemean,femalevariance),
    logit(femalemean) <-a+b*femalesparseness+cf*dominance_femalevalues,
    malesteepness ~ dbeta2(malemean,malevariance),
    logit(malemean) <-a+b*malesparseness+(cf+maleoffset)*dominance_malevalues,
    a~dnorm(0,1),
    b~dnorm(0,1),
    cf~dnorm(0,1),
    maleoffset~dnorm(0,1),
    femalevariance~dexp(10),
    malevariance~dexp(10)
  ) , data=mdata_sex_steepness_dominance , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)

# Check whether the male-offset is different from zero - if it is positive, that means that male hierarchies are steeper than female hierarchies when females are dominant; if it is negative, than male hierarchies are less steep than female hierarchies when females are dominant:
precis(m_sex_steepness_dominance)
# We can also check the cf value - if it is positive, that means that hierarchies of females are steeper when females win more fights and if the maleoffset is zero, it also means that male hierarchies are steeper when females win more fights; if it is negative, hierarchies are steeper when males are dominant.

# We can plot this - female values are black, male values are red
plot(data_dominance$steepness~data_dominance$perc_won_females,col=as.factor(data_dominance$sex))
# It looks like there is a confound here - we do not have hierarchy data for males from species where males always wins the fights. This is because the species where males always win fights, there is usually only a single male per group - that means we cannot calculate the hierarchy among the males (in our data, those species are gorillas, hamadryas baboons, and red howler monkeys, which are all polygynous). 


#PLOT C.10)
#predictions
dominance_seq <- seq(
  min(mdata_sex_steepness_dominance$dominance_femalevalues),
  max(mdata_sex_steepness_dominance$dominance_femalevalues),
  length.out = 100
)

#posteriors means
post <- extract.samples(m_sex_steepness_dominance)

pred_females <- inv_logit(mean(post$a) + mean(post$b)*0 + mean(post$cf)*dominance_seq)
pred_males <- inv_logit(mean(post$a) + mean(post$b)*0 + (mean(post$cf) + mean(post$maleoffset))*dominance_seq)

#df
plot_pred <- data.frame(
  dominance = rep(dominance_seq, 2),
  steepness = c(pred_females, pred_males),
  sex = rep(c("Females","Males"), each = length(dominance_seq))
)

#plot1 : steepness vs dominance femelle
C10p1 <- ggplot(plot_pred, aes(x = dominance, y = steepness, color = sex)) +
  geom_line(size = 1.5) +
  geom_point(
    data = data_dominance,
    aes(
      x = standardize(perc_won_females),
      y = steepness,
      color = sex
    ),
    alpha = 0.3
  ) +
  scale_color_manual(values=c("Females"="#443A83FF","Males"="gold3")) +
  xlab("Dominance femelle (standardisée)") +
  ylab("Steepness") +
  theme_minimal(base_size = 16) +
  theme(legend.position = "top")

C10p1


                           
#   12) Are hierarchies for females more likely to be based on signals whereas those in males more likely to be based on aggression?----

# summary:
data %>% group_by(sex,typeofbehaviour) %>% summarise(n())
# For both sexes, we have very few matrices that only contain data from aggressive interactions (10 of the 156)
# We therefore reclassify this into whether the dominance interactions involved any aggression at all (A + AD) or whether they are purely based on symbols (D)

data$binarybehaviour<-ifelse(data$typeofbehaviour=="D","D","A")

data %>%
  mutate(binarybehaviour = ifelse(typeofbehaviour == "D", "Signals (D)", "Agression (A+AD)")) %>%
  count(sex, binarybehaviour, sort = TRUE) %>%
  mutate(pourcentage = round(n / sum(n) * 100, 1))


mdata_sex_behaviour <- list(
  behaviour=as.numeric(as.factor(data$binarybehaviour))-1,
  sex=as.numeric(as.factor(data$sex))
)

m_sex_behaviour <- ulam(
  alist(
    behaviour ~ dbinom(1,p),
    logit(p) <-a[sex],
    a[sex]~dnorm(0,1)
  ) , data=mdata_sex_behaviour , chains=4 , cores=4 , cmdstan=T, messages=FALSE, refresh=0)

posterior_sex_behaviour<-extract.samples(m_sex_behaviour)
mean_prop_females<-inv_logit(posterior_sex_behaviour$a[,1])
mean_prop_males<-inv_logit(posterior_sex_behaviour$a[,2])
contrast_sex_behaviour<-inv_logit(posterior_sex_behaviour$a[,2])-inv_logit(posterior_sex_behaviour$a[,1])
precis(contrast_sex_behaviour)
results_sex_behaviour<-list(prop_onlydisplay_females=mean_prop_females,prop_onlydisplay_males=mean_prop_males,difference_prob_onlydisplay=contrast_sex_behaviour)

precis(results_sex_behaviour)


### Estimation of whether the phylogenetic component is captured by the number of interactions----
A<-vcv.phylo(mtree)

model_phy<-brm(
  steepness ~ 1 + (1|gr(species,cov=A)),
  data=data[data$sex=="females",],
  data2=list(A=A),
  family=Beta()
)
summary(model_phy)

model_phy_interactions<-brm(
  steepness ~ l_interactions+ (1|gr(species,cov=A)),
  data=data[data$sex=="females",],
  data2=list(A=A),
  family=Beta()
)
summary(model_phy_interactions)

# The estimate for the intercept of the phylogenetic covariance declines, but is still largely present when accounting for the number of interactions

save.image("~/Desktop/steepness_project/Analyse/RData_Analyse30-09.RData")
