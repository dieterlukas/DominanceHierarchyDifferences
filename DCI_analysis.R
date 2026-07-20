# This code is for the analyses of the dci in the the article by Spicher, Huchard, Lukas 
# The aim is to identify whether the dci of males and females differ among primates

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
library(ppcor)
library(ggcorrplot)
library(factoextra)


# Load the data
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
data <- read.csv2("Comparative_primate_steepness_values.csv")


data$dci <- as.numeric(data$dci)
data$h_index <- as.numeric(data$h_index)
data$steepness <- as.numeric(data$steepness)

mean(data$dci, na.rm = TRUE)
mean(data$h_index, na.rm = TRUE)
mean(data$steepness, na.rm = TRUE)

mean(data$dci[data$sex == "males"], na.rm = TRUE)
sd(data$dci[data$sex == "males"], na.rm = TRUE)

mean(data$dci[data$sex == "females"], na.rm = TRUE)
sd(data$dci[data$sex == "females"], na.rm = TRUE)


#We have a lot of dci = 1 values, which fails the convergence of the model, so prevent it from reaching the bound of 1 
data$dci2 <- data$dci
data$dci2[data$dci2 == 1] <- 0.9999

# Is DCI senstive to data sparseness? 
dat_list_dci_sparseness <- list(
  dci2 = as.numeric(data$dci2),  
  sparseness = as.numeric(data$sparseness)
)

m_dci_sparseness <- ulam(
  alist(
    dci2 ~ dbeta2(mean,variance),
    logit(mean) <-a + b*sparseness,
    a ~dnorm(0.5,1),   
    b ~dnorm(0,1),
    variance ~ dexp(10)
  ) , data=dat_list_dci_sparseness , chains=4 , cores=8 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)
# We check the results - we are interested in the effect measured in the factor b
precis(m_dci_sparseness)
plot(data$dci~data$sparseness)
#no effect of sparseness on the dci, as suggesting by Koenig et al., 2013 


# Is DCI senstive to the density of the matrice (i.e. number of interactions)? 
dat_list_dci_density <- list(
  dci2 = as.numeric(data$dci2),  
  density = log(data$numberofineractions)
)

m_dci_density <- ulam(
  alist(
    dci2 ~ dbeta2(mean,variance),
    logit(mean) <-a + b*density,
    a ~dnorm(0.5,1),   
    b ~dnorm(0,1),
    variance ~ dexp(10)
  ) , data=dat_list_dci_density , chains=4 , cores=8 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)
# We check the results - we are interested in the effect measured in the factor b
precis(m_dci_density)
plot(data$dci~log(data$numberofineractions))
#DCI decreases as the number of recorded interactions increases
#The 89% credible interval for b ranges from -0.24 to -0.04, so it is entirely negative. This suggests a reasonably credible effect


# Is DCI different in captivity vs in wild? 
dat_list_dci_captivity <- list(
  dci = as.numeric(data$dci2),
  captive = as.numeric(as.factor(data$group)),
  numberofindividuals= standardize(data$numberofindividuals)
)

# We need to provide priors, our expectation of what these values might be. For the means, we could expect that they are somewhere around 0.5
# For the variance, we expect this to be larger than than zero (so we use the dexp function) and larger than one
# For the relationship with captivity, we don't expect any effect a priori
m_dci_captivity <- ulam(
  alist(
    dci ~ dbeta2(mean,variance),
    logit(mean) <-a + b[captive]+c*numberofindividuals,
    a ~dnorm(0.5,1),   
    b[captive] ~dnorm(0,1),
    c ~dnorm(0,1),
    variance ~ dexp(10)
  ) , data=dat_list_dci_captivity , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)
# We check the results - we are interested in the effect measured in the factor b
posterior_captivity<-extract.samples(m_dci_captivity)
contrast_captivity<-inv_logit(posterior_captivity$b[,2])-inv_logit(posterior_captivity$b[,1])
precis(contrast_captivity)
#no effect




#PLOT figure 2 : is there a difference of dci between our main primate group? 

data$clade <- case_when(
  data$species %in% c("Pan_troglodytes", "Pan_paniscus", "Gorilla_beringei", "Gorilla_gorilla") ~ "Hominoidea",
  data$species %in% c("Cebus_capucinus", "Cebus_apella", "Saimiri_sciureus", "Alouatta_palliata") ~ "Platyrrhini",
  TRUE ~ "Cercopithecoidea"
)

#PLOT
# 1. Stats PAR CLADES (mâles+femelles)
clade_stats <- data %>%
  group_by(clade) %>%
  summarise(
    N_matrices = n(),
    N_species = n_distinct(species),
    
    # STEEPNESS
    median_dci = round(median(dci, na.rm=T), 3),
    mean_dci= round(mean(dci, na.rm=T), 3),
    ci_lower_steep = round(quantile(dci, 0.10, na.rm=T), 3),
    ci_upper_steep = round(quantile(dci, 0.90, na.rm=T), 3),
    )
print(clade_stats)

label_map <- c("Cercopithecoidea" = "Cercopithecoidea", 
               "Hominoidea" = "Hominoidea", 
               "Platyrrhini" = "Platyrrhini")

#plot steepness median
p_clades_dci_med <- ggplot() +
  geom_jitter(data = data,
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




#Are DCI and steepness correlated across species?

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


mdata_phylogeny_dci_steepness <- list(
  dci=data$dci2,
  steepness=data$steepness,
  species=as.integer(as.factor(data$species)),
  N_spp=length(unique(data$species))
)

Dmat<-cophenetic(mtree)
mdata_phylogeny_dci_steepness$Dmat<-Dmat/max(Dmat)
colnames(mdata_phylogeny_dci_steepness$Dmat)<-as.integer(as.factor(colnames(mdata_phylogeny_dci_steepness$Dmat)))
rownames(mdata_phylogeny_dci_steepness$Dmat)<-as.integer(as.factor(rownames(mdata_phylogeny_dci_steepness$Dmat)))

m_dci_steepness <- ulam(
  alist(
    dci ~ dbeta2(mean,variance),
    logit(mean) <-a+b[species]+c*steepness,
    a ~dnorm(0,1),
    c~dnorm(0,1),
    vector[N_spp]:b~multi_normal(0,SIGMA),
    matrix[N_spp,N_spp]:SIGMA <- cov_GPL2( Dmat , etasq, rhosq , 0.01 ),
    etasq ~ half_normal(1,0.25),
    rhosq ~ half_normal(3,0.25),
    variance~dexp(10)
  ) , data=mdata_phylogeny_dci_steepness , chains=4 , cores=8 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)

precis(m_dci_steepness)
plot(data$dci~data$steepness)
# DCI and steepness are positively correlated


#PLOT DCI - steepness - figure 3
data_plot_3b <- data.frame(
  steepness = mdata_phylogeny_dci_steepness$steepness,
  dci = mdata_phylogeny_dci_steepness$dci,
  species = factor(data$species)
)

post <- extract.samples(m_dci_steepness)

steepness_seq <- seq(min(data_plot_3b$steepness), max(data_plot_3b$steepness), length.out = 100)

preds <- sapply(seq_along(post$a), function(i) {
  plogis(post$a[i] + post$c[i] * steepness_seq)
})

pred_df <- data.frame(
  steepness = steepness_seq,
  mean = apply(preds, 1, mean),
  lower = apply(preds, 1, quantile, 0.055),
  upper = apply(preds, 1, quantile, 0.945)
)

#rawdata
data_plot_3b <- data.frame(
  dci = mdata_phylogeny_dci_steepness$dci,
  steepness = mdata_phylogeny_dci_steepness$steepness,
  species = factor(data$species)  
)

data_plot_3b$clade <- dplyr::case_when(
  data_plot_3b$species %in% c("Cebus_capucinus", "Saimiri_sciureus", "Alouatta_palliata") ~ "Platyrrhini",
  data_plot_3b$species %in% c("Pan_paniscus", "Pan_troglodytes", "Gorilla_gorilla", "Gorilla_beringei") ~ "Hominoidea",
  TRUE ~ "Cercopithecoidea"
)

ggplot(data_plot_3b, aes(x = steepness, y = dci, color = clade)) +
  
  geom_ribbon(data = pred_df,
              aes(x = steepness, ymin = lower, ymax = upper),
              fill = "#fcbba1", alpha = 0.20, inherit.aes = FALSE) +
  
  geom_line(data = pred_df,
            aes(x = steepness, y = mean),
            color = "#d73027", size = 1.2, inherit.aes = FALSE) +
  
  geom_point(size = 3, stroke = 0.4,
             position = position_jitter(width = 0, height = 0.01)) +
  
  scale_color_manual(values = c(
    Platyrrhini      = "#8da0cb",
    Hominoidea       = "#fc8d62",
    Cercopithecoidea = "#66c2a5"
  )) +
  
  labs(x = "Hierarchy steepness",
       y = "DCI",
       color = NULL) +
  
  theme_light(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    legend.position = c(0.02, 0.02),
    legend.justification = c(0, 0),
    legend.background = element_rect(fill = alpha("white", 0.9), color = NA),
    legend.text = element_text(size = 10)
  )




#Are DCI and linearity correlated across species?
mdata_phylogeny_dci_linearity <- list(
  dci=data$dci2,
  linearity=data$h_index,
  species=as.integer(as.factor(data$species)),
  N_spp=length(unique(data$species))
)

Dmat<-cophenetic(mtree)
mdata_phylogeny_dci_linearity$Dmat<-Dmat/max(Dmat)
colnames(mdata_phylogeny_dci_linearity$Dmat)<-as.integer(as.factor(colnames(mdata_phylogeny_dci_linearity$Dmat)))
rownames(mdata_phylogeny_dci_linearity$Dmat)<-as.integer(as.factor(rownames(mdata_phylogeny_dci_linearity$Dmat)))

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
  ) , data=mdata_phylogeny_dci_linearity , chains=4 , cores=4 , log_lik=TRUE , cmdstan=T, messages=FALSE, refresh=0)

precis(m_dci_linearity)
plot(data$dci~data$h_index)
# c crosses 0 : CI = [0.00; 1.25]


#PLOT DCI - linearity - figure 3
data_plot_3c <- data.frame(
  linearity = mdata_phylogeny_dci_linearity$linearity,
  dci = mdata_phylogeny_dci_linearity$dci,
  species = factor(data$species)
)

post <- extract.samples(m_dci_linearity)

linearity_seq <- seq(min(data_plot_3c$linearity), max(data_plot_3c$linearity), length.out = 100)

preds <- sapply(seq_along(post$a), function(i) {
  plogis(post$a[i] + post$c[i] * linearity_seq)
})

pred_df <- data.frame(
  linearity = linearity_seq,
  mean = apply(preds, 1, mean),
  lower = apply(preds, 1, quantile, 0.055),
  upper = apply(preds, 1, quantile, 0.945)
)

#rawdata
data_plot_3c <- data.frame(
  dci = mdata_phylogeny_dci_linearity$dci,
  linearity = mdata_phylogeny_dci_linearity$linearity,
  species = factor(data$species)  
)

data_plot_3c$clade <- dplyr::case_when(
  data_plot_3c$species %in% c("Cebus_capucinus", "Saimiri_sciureus", "Alouatta_palliata") ~ "Platyrrhini",
  data_plot_3c$species %in% c("Pan_paniscus", "Pan_troglodytes", "Gorilla_gorilla", "Gorilla_beringei") ~ "Hominoidea",
  TRUE ~ "Cercopithecoidea"
)

ggplot(data_plot_3c, aes(x = linearity, y = dci, color = clade)) +
  
  geom_ribbon(data = pred_df,
              aes(x = linearity, ymin = lower, ymax = upper),
              fill = "#fcbba1", alpha = 0.20, inherit.aes = FALSE) +
  
  geom_line(data = pred_df,
            aes(x = linearity, y = mean),
            color = "#d73027", size = 1.2, inherit.aes = FALSE) +
  
  geom_point(size = 3, stroke = 0.4,
             position = position_jitter(width = 0, height = 0.01)) +
  
  scale_color_manual(values = c(
    Platyrrhini      = "#8da0cb",
    Hominoidea       = "#fc8d62",
    Cercopithecoidea = "#66c2a5"
  )) +
  
  labs(x = "Hierarchy linearity",
       y = "DCI",
       color = NULL) +
  
  theme_light(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    legend.position = c(0.02, 0.02),
    legend.justification = c(0, 0),
    legend.background = element_rect(fill = alpha("white", 0.9), color = NA),
    legend.text = element_text(size = 10)
  )



df_corr <- data.frame(
  steepness = mdata_phylogeny_dci_steepness$steepness,
  linearity = mdata_phylogeny_dci_linearity$linearity,
  dci = mdata_phylogeny_dci_steepness$dci
)

# Pearson’s correlations
cor_mat <- cor(df_corr, use = "pairwise.complete.obs", method = "pearson")
round(cor_mat, 2)

ggcorrplot(cor_mat,
           type = "full",
           method = "square",
           hc.order = TRUE,
           lab = TRUE,
           lab_size = 5,
           colors = c("#2166ac", "white", "#b2182b"),
           outline.col = "grey90",
           title = "Correlation matrix",
           ggtheme = theme_minimal(base_size = 14)) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
#The correlation matrix shows a strong positive association between DCI and steepness, whereas DCI and linearity are only weakly correlated. 
#Steepness and linearity are moderately correlated, suggesting partial overlap but not redundancy.




df_pca <- data.frame(
  steepness = mdata_phylogeny_dci_steepness$steepness,
  linearity = mdata_phylogeny_dci_linearity$linearity,
  dci = mdata_phylogeny_dci_steepness$dci
) %>%
  drop_na()

# PCA
res_pca <- prcomp(df_pca, scale. = FALSE)

# PCA results
summary(res_pca)
res_pca$rotation
res_pca$x
#A PCA on steepness, linearity, and DCI revealed a strong shared component, with the first principal component explaining 59.4% of the variance and loading positively on all three variables.
#However, a second component explained an additional 36.8% and mainly contrasted linearity against DCI and steepness, indicating that the three metrics are related but not fully redundant.

res_pca$rotation

fviz_eig(res_pca,
         addlabels = TRUE,
         barfill = "#4C78A8",
         barcolor = "white",
         linecolor = "#1f4e79") +
  theme_minimal(base_size = 14) +
  labs(title = "Variance explained by PCA axes")

# 2. Correlation Circle
fviz_pca_var(res_pca,
             col.var = "contrib",
             gradient.cols = c("#9ecae1", "#f7f7f7", "#d73027"),
             repel = TRUE) +
  theme_minimal(base_size = 14) +
  labs(title = "PCA")

# 3. Biplot 
fviz_pca_biplot(res_pca,
                label = "var",
                col.var = "#d73027",
                col.ind = "grey70",
                alpha.ind = 0.6,
                repel = TRUE) +
  theme_minimal(base_size = 14) +
  labs(title = "PCA biplot")

#DCI and steepness appear to reflect the same aspect of the hierarchy, whilst linearity shares this general structure but also provides slightly different information. 
#This is consistent with previous correlations: a strong relationship between DCI and steepness, and a weaker relationship between DCI and linearity.



#Are DCI and steepness still linked if we remove the effect of linearity?



df_pcor <- data.frame(
  steepness = mdata_phylogeny_dci_steepness$steepness,
  linearity  = mdata_phylogeny_dci_linearity$linearity,
  dci        = mdata_phylogeny_dci_steepness$dci
) %>%
  drop_na()

pcor_res <- pcor(df_pcor, method = "pearson")
pcor_res
#Partial correlations showed that the association between DCI and steepness remained strong after controlling for linearity, whereas the residual association between DCI and linearity was weak and negative. 
#This suggests that steepness captures the main shared component linking the dominance metrics






# Is DCI in males different from that in females (w. or without accounting for phylogenetic relatedness)?
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

#2nd plot
cols_sex <- c("males"="gold3", "females"="#443A83FF")

ggplot(data, aes(x = sex, y = dci, color = sex, fill = sex)) +
  geom_violin(alpha = 0.2, width = 0.3, color = NA) +
  stat_summary(fun = median, geom = "point", size = 4, color = "black") +
  geom_jitter(width = 0.15, size = 3, alpha = 0.6, shape = 21, stroke = 0.5) +
  scale_color_manual(values = cols_sex) +
  scale_fill_manual(values = cols_sex) +
  labs(x = "Sex", y = "DCI", title = NULL) +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    legend.position = c(0, 1),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = scales::alpha("white", 0.8), color = NA)
  )

#3rd plot
height_points <- c("males" = 6, "females" = 7.5)

steep_descri <- ggplot(data, aes(x = dci, color = sex, fill = sex)) +
  geom_density(size = 2, alpha = 0.3) +
  geom_jitter(data = data,
              aes(y = height_points[sex]),
              width = 0, size = 3, alpha = 0.6, shape = 21, fill = "white", stroke = 1.5) +
  scale_color_manual(values = cols_sex) +
  scale_fill_manual(values = cols_sex) +
  labs(x = "DCI", y = "Density") +
  xlim(0,1.04) +
  ylim(0,8) + 
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face="bold", hjust=0.5),
    axis.title = element_text(face="bold", size=16),
    axis.text = element_text(size=14),
    legend.position="top"
  )

steep_descri



data$sexspecies<-paste(data$sex,data$species,sep="_")

dat_list_dci_species <- list(
  dci = data$dci2,  
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


##############################################################################
# Phylogenetic analyses

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



mdata_phylogeny_both <- list(
  dci_females=dcidata[dcidata$sex=="females",]$dci2,
  species_females=as.integer(as.factor((dcidata[dcidata$sex=="females",]$species))),
  N_spp_females=length(unique(dcidata[dcidata$sex=="females",]$species)),
  dci_males=dcidata[dcidata$sex=="males",]$dci2,
  species_males=as.integer(as.factor((dcidata[dcidata$sex=="males",]$species))),
  N_spp_males=length(unique(dcidata[dcidata$sex=="males",]$species))
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















