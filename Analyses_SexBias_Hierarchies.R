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

# Load the data
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
data <- read.csv("Comparative_primate_steepness_values.csv")

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
data[data$studyperiod_month=="unc",]$studyperiod_month<-NA
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



# 1b) steepness sensitive to number of interactions
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


# 1c) linearity sensitive to sparseness
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


# 1d) linearity sensitive to number of interactions
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

# 2a) steepness sensitive to the nb of individuals
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


# 2b) linearity sensitive to the nb of individuals
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
### DIETER there is an effect here!!


# 3a) steepness different in captivity than in the wild
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


# 3b) linearity different in captivity than in the wild
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

#TREE
p <- ggtree(mtree) + theme_tree2()

data_N <- data %>%
  group_by(species) %>%
  summarise(Nmatrices = n(), .groups = "drop") %>%  # n() = nombre de lignes par espèce
  mutate(
    species_label = str_replace_all(species, "_", " "),
    species_label = str_to_sentence(species_label)
  )
p2 <- p %<+% data_N

p2 +
  geom_hilight(node = 70, fill = "#8da0cb", alpha = 0.2) +
  geom_hilight(node = 41, fill = "#66c2a5", alpha = 0.2) +
  geom_hilight(node = 67, fill = "#fc8d62", alpha = 0.2) +
  
  geom_tippoint(aes(size = Nmatrices), shape = 21, fill = "steelblue", alpha = 0.7) +
  
  geom_tiplab(aes(label = species_label),
              size = 3.5, offset = 4, align = TRUE,
              fontface = "italic") +
  
  geom_cladelabel(node = 70, label = "Platyrrhini",
                  align = TRUE,
                  offset = 25,
                  offset.text = 3,
                  fontsize = 4,
                  barsize = 1) +
  
  geom_cladelabel(node = 41, label = "Cercopithecoidea",
                  align = TRUE,
                  offset = 25,
                  offset.text = 3,
                  fontsize = 4,
                  barsize = 1) +
  
  geom_cladelabel(node = 67, label = "Hominoidea",
                  align = TRUE,
                  offset = 25,
                  offset.text = 3,
                  fontsize = 4,
                  barsize = 1) +
  scale_size_area(
    max_size = 6,
    name = "N matrices",
    breaks = round(seq(min(data_N$Nmatrices),
                       max(data_N$Nmatrices),
                       length.out = 4))
  ) +
  
  coord_cartesian(xlim = c(0, 95), clip = "off")





### Analyses for section B - dominance types

#   4) Are hierarchy steepness and linearity correlated across species  / accounting for phylogenetic relatedness and sparseness (which affects both)

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

###############################################################################################
# PLOT for B.4) 

##1st plot -  classic

#posterior
post <- extract.samples(m_steepness_steepness_linearity)

#smooth curve
linearity_seq <- seq(min(plot_4_data$linearity), max(plot_4_data$linearity), length.out = 100)

#predictions
preds <- sapply(1:length(post$a), function(i) {
  plogis(post$a[i] + post$d[i] * linearity_seq)  # steepness sur logit scale
})

#mean and CI
pred_df <- data.frame(
  linearity = linearity_seq,
  mean = rowMeans(preds),
  lower = apply(preds, 1, quantile, 0.055),
  upper = apply(preds, 1, quantile, 0.945)
)

#rawdata
plot_4_data <- data.frame(
  steepness = mdata_phylogeny_steepness_linearity$steepness,
  linearity = mdata_phylogeny_steepness_linearity$linearity,
  species = factor(data$species)  # if we want to put some colours
)



### plot classic
ggplot(plot_4_data, aes(x = linearity, y = steepness)) +
  geom_point(size = 3, color = "steelblue") +
  geom_ribbon(data = pred_df, aes(x = linearity, ymin = lower, ymax = upper),
              fill = "red", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = pred_df, aes(x = linearity, y = mean),
            color = "red", size = 1.2, inherit.aes = FALSE) +
  labs(
    x = "Hierarchy Linearity",
    y = "Hierarchy Steepness") +
  theme_bw(base_size = 14) +
  theme(axis.title = element_text(face = "bold"))


### color = sex
plot_4_data_sex <- data.frame(
  steepness = mdata_phylogeny_steepness_linearity$steepness,
  linearity = mdata_phylogeny_steepness_linearity$linearity,
  species = factor(data$species),
  sex = factor(data$sex)
)

ggplot(plot_4_data_sex, aes(x = linearity, y = steepness, color = sex)) +
  
  geom_ribbon(data = pred_df,
              aes(x = linearity, ymin = lower, ymax = upper),
              fill = "#fcbba1",
              alpha = 0.3,
              inherit.aes = FALSE) +
  
  geom_line(data = pred_df,
            aes(x = linearity, y = mean),
            color = "#d73027",
            size = 1.2,
            inherit.aes = FALSE) +
  
  geom_point(size = 3) +
  
  scale_color_manual(values = c(
    females = "#443A83FF",
    males = "#FDE725FF"
  )) +
  
  labs(
    x = "Hierarchy Linearity",
    y = "Hierarchy Steepness",
    color = NULL
  ) +
  
  theme_bw(base_size = 14) +
  
  theme(
    axis.title = element_text(face = "bold"),
    legend.position = c(0.005, 0.05),
    legend.justification = c(0, 0.3),
    legend.background = element_rect(fill = alpha("white", 0.8), color = NA)
  )


### color = sex and shape = clade
plot_4_data_sex$clade <- dplyr::case_when(
  plot_4_data_sex$species %in% c("Cebus_capucinus", "Saimiri_sciureus", "Alouatta_palliata") ~ "Platyrrhini",
  plot_4_data_sex$species %in% c("Pan_paniscus", "Pan_troglodytes", "Gorilla_gorilla", "Gorilla_beringei") ~ "Hominoidea",
  TRUE ~ "Cercopithecoidea"
)

ggplot(plot_4_data_sex, aes(x = linearity, y = steepness)) +
  geom_ribbon(
    data = pred_df,
    aes(x = linearity, ymin = lower, ymax = upper),
    fill  = "#fcbba1",
    alpha = 0.20,
    inherit.aes = FALSE
  ) +
  geom_line(
    data = pred_df,
    aes(x = linearity, y = mean),
    color = "#d73027",
    size  = 1.2,
    inherit.aes = FALSE
  ) +
  geom_point(
    aes(color = sex, shape = clade),
    size = 3,
    stroke = 0.4,
    position = position_jitter(width = 0.01, height = 0)
  ) +
  scale_color_manual(values = c(
    females = "#443A83FF",
    males   = "#FDE725FF"
  )) +
  scale_shape_manual(values = c(
    Platyrrhini      = 16,
    Hominoidea       = 17,
    Cercopithecoidea = 15
  )) +
  labs(
    x = "Hierarchy Linearity",
    y = "Hierarchy Steepness",
    color = NULL,
    shape = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major   = element_line(color = "grey85"),
    panel.grid.minor   = element_blank(),
    axis.title         = element_text(face = "bold"),
    legend.background  = element_rect(fill = alpha("white", 0.9), color = NA),
    legend.position    = c(0.01, 0.01),
    legend.justification = c(0, 0),
    legend.text        = element_text(size = 9),
    legend.title       = element_text(size = 9),
    legend.box          = "horizontal",
    legend.box.spacing  = grid::unit(0.3, "cm"),,
    
  )

### color = clades
ggplot(plot_4_data_sex, aes(x = linearity, y = steepness, color = clade)) +
  
  geom_ribbon(data = pred_df,
              aes(x = linearity, ymin = lower, ymax = upper),
              fill = "#fcbba1", alpha = 0.20, inherit.aes = FALSE) +
  
  geom_line(data = pred_df,
            aes(x = linearity, y = mean),
            color = "#d73027", size = 1.2, inherit.aes = FALSE) +
  
  geom_point(size = 3, stroke = 0.4,
             position = position_jitter(width = 0.01, height = 0)) +
  
  scale_color_manual(values = c(
    Platyrrhini      = "#8da0cb", 
    Hominoidea       = "#fc8d62",    
    Cercopithecoidea = "#66c2a5"     
  )) +
  
  labs(x = "Hierarchy Linearity", 
       y = "Hierarchy Steepness",
       color = NULL) +
  
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major   = element_line(color = "grey85"),
    panel.grid.minor   = element_blank(),
    axis.title         = element_text(face = "bold"),
    legend.position    = c(0.02, 0.05),
    legend.justification = c(0, 0),
    legend.background  = element_rect(fill = alpha("white", 0.9), color = NA),
    legend.text        = element_text(size = 10)
  )


#   5) Are hierarchies steeper and more linear when they are based on signals rather than aggression?
data$typeofbehaviour[is.na(data$typeofbehaviour)] <- "AD"
#if we don't want to include NA : data <- data[!is.na(data$typeofbehaviour), ]

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

#PLOTs for section B.5) 
# We need that typeofbehaviour is a factor
data$typeofbehaviour <- factor(data$typeofbehaviour, levels = c("A","D","AD"))

ggplot(data, aes(x = typeofbehaviour, y = steepness, fill = typeofbehaviour)) +
  geom_violin(alpha = 0.3, color = NA) +          # distribution générale
  geom_boxplot(width = 0.2, outlier.shape = 21, outlier.fill = "white") + # médiane + quartiles
  scale_fill_manual(values = c("A" = "#1F77B4", "D" = "#2E8B57", "AD" = "#FF7F0E")) +
  labs(
    x = "Type of Interaction",
    y = "Hierarchy Steepness") +
  theme_bw(base_size = 16) +
  theme(
    axis.text = element_text(size = 12),
    legend.position = "none"
  )



########################################################################################################
#PLOTs for section B.5) 

palette_complete <- wes_palette("GrandBudapest1")
print(palette_complete)          
print(tail(palette_complete))  

p_steepness <- ggplot() +
  geom_jitter(data = data,
              aes(x = typeofbehaviour, y = steepness, color = typeofbehaviour),
              width = 0.25, size = 5, alpha = 0.6, stroke = 0) +
  geom_errorbar(data = cred_steepness,
                aes(x = typeofbehaviour, ymin = pmax(0, lower), ymax = pmin(1, upper)),
                width = 0.25, size = 0.5, color = "black") +
  geom_point(data = cred_steepness,
             aes(x = typeofbehaviour, y = median),
             size = 5, shape = 21, fill = "white", color = "black", stroke = 1.2) +
  scale_color_manual(values = c("#F1BB7B", "#5B1A18", "#FD6467")) +  
  scale_x_discrete(labels = label_map) +
  labs(x = NULL, y = "Hierarchy Steepness") +
  ylim(0, 1) +
  theme_minimal(base_size = 16) +
  theme(legend.position = "none")


p_linearity <- ggplot() +
  geom_jitter(data = data,
              aes(x = typeofbehaviour, y = steepness, color = typeofbehaviour),
              width = 0.25, size = 5, alpha = 0.6, stroke = 0) +  
  geom_errorbar(data = cred_steepness,
                aes(x = typeofbehaviour, ymin = pmax(0, lower), ymax = pmin(1, upper)),
                width = 0.25,    # plus étroites
                size = 0.5,      # plus fines
                color = "black")+
  geom_point(data = cred_linearity,
             aes(x = typeofbehaviour, y = median),
             size = 5, shape = 21, fill = "white", color = "black", stroke = 1.2) +
  scale_color_manual(values = c("#F1BB7B", "#5B1A18", "#FD6467")) +  
  scale_x_discrete(labels = label_map) +
  labs(x = "Type of Interaction", y = "Hierarchy Linearity") +
  ylim(0, 1) +
  theme_minimal(base_size = 16) +
  theme(
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

# merge them
x.grob <- textGrob("Type of Interaction", gp = gpar(fontsize = 16, fontface = "bold"))
grid.arrange(p_steepness + xlab(NULL), p_linearity + xlab(NULL),
             ncol = 2, bottom = x.grob)



### Analyses for section C - sex-differences 
#   6) Is hierarchy steepness in males different from that in females (w. or without accounting for phylogenetic relatedness)?

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
cols_sex <- c("males"="#FDE725FF", "females"="#443A83FF")

ggplot(data, aes(x = sex, y = steepness, color = sex, fill = sex)) +
  geom_violin(alpha = 0.2, width = 0.3, color = NA) +
  stat_summary(fun = median, geom = "point", size = 4, color = "black") +
  geom_jitter(width = 0.15, size = 3, alpha = 0.6, shape = 21, stroke = 0.5) +
  scale_color_manual(values = cols_sex) +
  scale_fill_manual(values = cols_sex) +
  labs(x = "Sex", y = "Hierarchy Steepness", title = "Steepness by Sex") +
  ylim(0, 1) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none"
  )


#3rd plot
height_points <- c("males" = 4.5, "females" = 5.5)

ggplot(data, aes(x = steepness, color = sex, fill = sex)) +
  geom_density(size = 2, alpha = 0.3) +
  geom_jitter(data = data,
              aes(y = height_points[sex]),
              width = 0, size = 3, alpha = 0.6, shape = 21, fill = "white", stroke = 1.5) +
  scale_color_manual(values = cols_sex) +
  scale_fill_manual(values = cols_sex) +
  labs(x = "Hierarchy Steepness", y = "Density") +
  xlim(0,1) +
  ylim(0,6) + 
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face="bold", hjust=0.5),
    axis.title = element_text(face="bold", size=16),
    axis.text = element_text(size=14),
    legend.position="top"
  )

  

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



##############################################################################
# Phylogenetic analyses

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


###############################################################################################################                                                                                       
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
  scale_fill_manual(values = c("Males" = "#FDE725FF", "Females" = "#443A83FF")) +
  scale_color_manual(values = c("Males" = "#FDE725FF", "Females" = "#443A83FF")) +
  xlim(0, 1) +
  xlab("Hierarchy steepness") +
  ylab("Posterior density") +
  theme_bw(base_size = 16) +
  theme(
    legend.title = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold")
  ) 


                                                                                       
### Checking sex differences within species

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

################################################################################################################################################
#BONUS plot
plot_df <- data.frame(
  diff = dat_list_sexdifference$steepness_difference,
  type = factor(dat_list_sexdifference$within,
                levels = c(0,1),
                labels = c("Between species (random pairs)",
                           "Within species (same species)"))
)

post <- extract.samples(m_sexdifference)

b_post <- data.frame(b = post$b)

P1 <- ggplot(plot_df, aes(x = diff, fill = type, color = type)) +
  geom_density(alpha = 0.3, size = 1.2) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1) +
  scale_fill_manual(values = c("#3B528BFF", "#FDE725FF")) +
  scale_color_manual(values = c("#3B528BFF", "#FDE725FF")) +
  xlab("Sex difference in steepness (male − female)") +
  ylab("Density") +
  ggtitle("A. Sex differences in steepness") +
  theme_minimal(base_size = 15) +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    plot.title = element_text(face = "bold")
  )

P2 <- ggplot(b_post, aes(x = b)) +
  geom_density(fill = "#21908CFF", alpha = 0.4, linewidth = 1.2) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1) +
  xlab("Effect of species identity (b)") +
  ylab("Posterior density") +
  ggtitle("B. Estimate of species effect") +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold")
  )

grid.arrange(P1, P2, ncol = 2)



# 7) If steepness differs in males and females, can it be linked to the fact that female hierarchies often include more individuals?

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


                                                                                       
### PLOT for C.7)
individuals <- standardize(data$numberofindividuals)

dat_list_steepness_individuals_n <- list(
  steepness = c(female_steepness, male_steepness),
  sex = c(rep(1, N_female_observations), rep(2, N_male_observations)),
  individuals = c(individuals[data$sex=="females"], individuals[data$sex=="males"])
)

individuals_seq <- seq(
  min(dat_list_steepness_individuals_n$individuals, na.rm = TRUE),
  max(dat_list_steepness_individuals_n$individuals, na.rm = TRUE),
  length.out = 100
)


#predicted steepness
pred_females_ind <- inv_logit(mean(post_steepness_individuals$a[,1]) + mean(post_steepness_individuals$b) * individuals_seq)
pred_males_ind   <- inv_logit(mean(post_steepness_individuals$a[,2]) + mean(post_steepness_individuals$b) * individuals_seq)


#long
plot_pred_ind <- data.frame(
  individuals = rep(individuals_seq, 2),
  steepness = c(pred_females_ind, pred_males_ind),
  sex = rep(c("Females","Males"), each = length(individuals_seq))
)

#plot steepness ~ nombre d'individus
C7p1 <- ggplot(plot_pred_ind, aes(x=individuals, y=steepness, color=sex)) +
  geom_line(size=1.5) +
  geom_point(data=data, 
             aes(x=standardize(numberofindividuals), y=steepness, color=sex), 
             alpha=0.3) +
  scale_color_manual(values=c("Females"="#443A83FF","Males"="#FDE725FF")) +
  xlab("Standardized number of individuals") +
  ylab("Steepness") +
  theme_minimal(base_size=16) +
  theme(legend.position="top")



individuals <- standardize(data$numberofindividuals)

#mean
individuals_mean <- mean(dat_list_steepness_individuals_n$individuals)

#predicted steepness
pred_females <- inv_logit(post_steepness_individuals$a[,1] + post_steepness_individuals$b * individuals_mean)
pred_males   <- inv_logit(post_steepness_individuals$a[,2] + post_steepness_individuals$b * individuals_mean)

#data frame for ggplot
plot_post_model_ind <- data.frame(
  steepness = c(pred_females, pred_males),
  sex = rep(c("Females","Males"), each = length(pred_females))
)

#density plot
C7p2<-ggplot(plot_post_model_ind, aes(x = steepness, color = sex, fill = sex)) +
  geom_density(alpha = 0.3, size = 1) +
  scale_color_manual(values = c("Females"="#443A83FF", "Males"="#FDE725FF")) +
  scale_fill_manual(values = c("Females"="#443A83FF", "Males"="#FDE725FF")) +
  xlab("Steepness (posterior predictive)") +
  ylab("Density") +
  theme_minimal(base_size = 16) +
  theme(legend.position = "top")


#merge them
grid.arrange(C7p1, C7p2, ncol = 2)

                                                                                       

### 8) Is hierarchy linearity in males different from that in females (w. or without accounting for phylogenetic relatedness)?
data$h_index<-data$h_index-0.001

mean_males <- mean(as.numeric(data$h_index[data$sex == "males"]), na.rm = TRUE)
mean_females <- mean(as.numeric(data$h_index[data$sex == "females"]), na.rm = TRUE)
mean_difference <- mean_males - mean_females
# In the raw data, the average h_index values are very close: males 0.81, females 0.81, difference 0.00

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
cols_sex <- c("males"="#FDE725FF", "females"="#443A83FF")

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

ggplot(data, aes(x = steepness, color = sex, fill = sex)) +
  geom_density(size = 2, alpha = 0.3) +
  geom_jitter(data = data,
              aes(y = height_points[sex]),
              width = 0, size = 3, alpha = 0.6, shape = 21, fill = "white", stroke = 1.5) +
  scale_color_manual(values = cols_sex) +
  scale_fill_manual(values = cols_sex) +
  labs(x = "Hierarchy Steepness", y = "Density") +
  xlim(0,1) +
  ylim(0,6) + 
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face="bold", hjust=0.5),
    axis.title = element_text(face="bold", size=16),
    axis.text = element_text(size=14),
    legend.position="top"
  )                                                                         

                                                                                       

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



##############################################################################
# Phylogenetic analyses

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

#PLOT for C.8)
##############################################################################################                                                                                       
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
  scale_fill_manual(values = c("Males" = "#FDE725FF", "Females" = "#443A83FF")) +
  scale_color_manual(values = c("Males" = "#FDE725FF", "Females" = "#443A83FF")) +
  xlim(0, 1) +
  xlab("Hierarchy linearity") +
  ylab("Posterior density") +
  theme_bw(base_size = 16) +
  theme(
    legend.title = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold")
  ) 


#   9) Is hierarchy linearity linked to hierarchy steepness in the same way in males and in females?



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
  scale_color_manual(values=c("females"="#443A83FF", "males"="#FDE725FF")) +
  xlab("Linearity (h_index)") +
  ylab("Steepness") +
  theme_minimal(base_size = 16) +
  theme(legend.position="top")


#Combined plot
#rawdata
data_plot <- data.frame(
  steepness = data$steepness,
  h_index = data$h_index,
  sex = data$sex
)

#predictions
h_index_seq <- seq(min(data$h_index, na.rm=TRUE),
                   max(data$h_index, na.rm=TRUE),
                   length.out=100)
post <- extract.samples(m_sex_steepness_linearity)

#predicted steepness
pred_females_mat <- sapply(h_index_seq, function(h) inv_logit(post$af + post$b*0 + post$cf*h))
pred_males_mat   <- sapply(h_index_seq, function(h) inv_logit(post$am + post$b*0 + (post$cf + post$maleoffset)*h))

#median and CI
pred_females_df <- data.frame(
  h_index = h_index_seq,
  median = apply(pred_females_mat, 2, median),
  lower  = apply(pred_females_mat, 2, quantile, probs=0.055),
  upper  = apply(pred_females_mat, 2, quantile, probs=0.945),
  sex = "Females"
)

pred_males_df <- data.frame(
  h_index = h_index_seq,
  median = apply(pred_males_mat, 2, median),
  lower  = apply(pred_males_mat, 2, quantile, probs=0.055),
  upper  = apply(pred_males_mat, 2, quantile, probs=0.945),
  sex = "Males"
)

pred_plot <- rbind(pred_females_df, pred_males_df)

#plot steepness ~ linearity
C9_plot1 <- ggplot() +
  geom_point(data=data, aes(x=h_index, y=steepness, color=sex), alpha=0.5, size=3) +
  geom_ribbon(data=pred_plot, aes(x=h_index, ymin=lower, ymax=upper, fill=sex), alpha=0.2) +
  geom_line(data=pred_plot, aes(x=h_index, y=median, color=sex), size=1.5) +
  scale_color_manual(values=c("Females"="#443A83FF", "Males"="#FDE725FF")) +
  scale_fill_manual(values=c("Females"="#443A83FF", "Males"="#FDE725FF")) +
  xlab("Linearity (h_index)") +
  ylab("Steepness") +
  theme_minimal(base_size=16) +
  theme(legend.position="top")



# density plot
slopes_df <- data.frame(
  slope = c(post$cf, post$cf + post$maleoffset),
  sex = rep(c("Females", "Males"), each = nrow(post))
)

C9_plot2 <- ggplot(slopes_df, aes(x=slope, fill=sex, color=sex)) +
  geom_density(alpha=0.3, size=1) +
  scale_color_manual(values=c("Females"="#443A83FF", "Males"="#FDE725FF")) +
  scale_fill_manual(values=c("Females"="#443A83FF", "Males"="#FDE725FF")) +
  xlab("Posterior slope (linearity effect on steepness)") +
  ylab("Density") +
  theme_minimal(base_size=16) +
  theme(legend.position="top")

grid.arrange(C9_plot1, C9_plot2, ncol=2)
                           
#####################################################################
#PLOT steepness and linearity

# 1. Graphs without legend
steep_no_leg <- steep_descri + theme(legend.position = "none")
linearity_no_leg <- linearity_descri + theme(legend.position = "none")

# 2. extract legend
leg <- get_legend(steep_descri)

# 3. merge them
grid.arrange(
  leg,                                    
  arrangeGrob(steep_no_leg, linearity_no_leg, ncol = 2), 
  ncol = 1, 
  heights = c(0.15, 1)                    
)
#####################################################################

                           
#   10) Are hierarchies steeper for the sex that wins more fights - that is, is the proportion of intersexual fights that females win negatively related to the hierarchy steepness in males, and positively to the hierarchy steepness in females?

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
  scale_color_manual(values=c("Females"="#443A83FF","Males"="#FDE725FF")) +
  xlab("Dominance femelle (standardisée)") +
  ylab("Steepness") +
  theme_minimal(base_size = 16) +
  theme(legend.position = "top")

C10p1


                           
#   11) Are hierarchies for females more likely to be based on signals whereas those in males more likely to be based on aggression?

# summary:
data %>% group_by(sex,typeofbehaviour) %>% summarise(n())
# For both sexes, we have very few matrices that only contain data from aggressive interactions (10 of the 156)
# We therefore reclassify this into whether the dominance interactions involved any aggression at all (A + AD) or whether they are purely based on symbols (D)

data$binarybehaviour<-ifelse(data$typeofbehaviour=="D","D","A")

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




### Estimation of whether the phylogenetic component is captured by the number of interactions
library(brms)
library(ape)
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


### Estimation of the amount of within compared to between species variance in steepness
dat_list_variation<-list(
  steepness=data[data$sex=="females",]$steepness,
  species=as.numeric(as.factor(data[data$sex=="females",]$species))
)

model <- ulam(
  alist(
    steepness ~ dbeta(mu, phi),
    logit(mu) <- a_group[species],
    a_group[species] ~ dnorm(a, sigma_group),
    a ~ dnorm(0, 1),
    sigma_group ~ dcauchy(0, 1),
    phi ~ dexp(1)
  ),
  data = dat_list_variation,
  chains = 4,
  cores = 4
)


# Extract posterior samples
post <- extract.samples(model)

post$mu<-inv_logit(post$a)
post$within<-post$mu*(1-post$mu)/post$phi
post$between<-post$sigma_group*post$sigma_group

# Compare variances
var_ratio <- post$between / post$within
precis(list(ratio=var_ratio))


dat_list_variation<-list(
  steepness=data[data$sex=="males",]$steepness,
  species=as.numeric(as.factor(data[data$sex=="males",]$species))
)

model <- ulam(
  alist(
    steepness ~ dbeta(mu, phi),
    logit(mu) <- a_group[species],
    a_group[species] ~ dnorm(a, sigma_group),
    a ~ dnorm(0, 1),
    sigma_group ~ dcauchy(0, 1),
    phi ~ dexp(1)
  ),
  data = dat_list_variation,
  chains = 4,
  cores = 4
)


# Extract posterior samples
post <- extract.samples(model)

post$mu<-inv_logit(post$a)
post$within<-post$mu*(1-post$mu)/post$phi
post$between<-post$sigma_group*post$sigma_group

# Compare variances
var_ratio <- post$between / post$within
precis(list(ratio=var_ratio))

op<-par()
plot.new()
par(mar = c(14.1, 4.1, 4.1, 4.1), # change the margins
    lwd = 2, # increase the line thickness
    cex.axis = 1.2 # increase default axis label size
)

plot(data[data$sex=="females",]$steepness~factor(data[data$sex=="females",]$species,levels=mtree$tip.label),las=2,xlab="",ylab="")
mtext("Hierarchy steepness in females",side=3,cex=2)

plot.new()
plot(data[data$sex=="males",]$steepness~factor(data[data$sex=="males",]$species,levels=mtree$tip.label),las=2,xlab="",ylab="")
mtext("Hierarchy steepness in males",side=3,cex=2)
