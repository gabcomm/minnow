#code written by Dr. Michael D. Burns (mdburns@ucdavis.edu)

require(phytools)

tree<-read.nexus("tree.nexus")
dat<-read.csv("data_reduced.csv")

dd<-name.check(tree, dat, data.names=dat$Species)
tree.final<-drop.tip(tree,dd$tree_not_data)

trop.func<-dat$Trophic
names(trop.func)<-dat$Species

func_simmaps <- make.simmap(tree = tree.final, x=trop.func, Q="empirical", pi = "estimated", nsim = 1000)

simmap.function.desc<-describe.simmap(func_simmaps)

cols<-setNames(palette()[1:6],mapped.states(func_simmaps[1]))

pdf(file="FigureXX_simmap_cypriniformes_fan.pdf", height=8, width=8)
plot(simmap.function.desc, colors=cols, fsize=0.6, cex=c(0.3,0.2), type="fan")
add.simmap.legend(colors=cols, prompt=FALSE,x=3, y=90, fsize=0.7)
dev.off()


##################################################
#         OUWie BM OU1
#################################################
library("tidyverse")

###################  Set up data frame for results #############################

best_model_data_2 <- tibble(.rows = 10000) # 5 traits * 2 models * 1000 simmaps
best_model_data_2$trait <- rep(NA,10000)
best_model_data_2$simmap_count <- rep(NA,10000)
best_model_data_2$model <- rep(NA,10000)
best_model_data_2$loglik <- rep(NA,10000)
best_model_data_2$aicc <- rep(NA,10000)
best_model_data_2$eigval <- rep(NA, 10000)
best_model_data_2$sigma.sq <- rep(NA, 10000)
best_model_data_2$alpha <- rep(NA,10000)
best_model_data_2$theta <- rep(NA,10000)
best_model_data_2$theta_se <- rep(NA,10000) #standard error
best_model_data_2$saddle <- rep(NA,10000)

#########################  Run OUwie  ##########################################

library(OUwie)

# define function
OUwie.model <- function(model, phy, data) {
  print(paste("Now starting model: ", model))
  return(OUwie(phy, data, model, simmap.tree=T,algorithm="invert", diagn=T)) 
}

models <- c("BM1","OU1")

trait_vector <- c("IL_log","IL_resid","ABV_log","ABV_resid", "Size")

j = 1 #trait
i = 1 #simmap
r = 1 
q = 1 #model

# subset data frame for each trait
for (j in 1:5){
  print(paste("Now starting trait: ", trait_vector[j]))
  OUwie_data <- data.frame(dat[ ,c(1,2,j+2)]) 
  i = 1
  
  # extract simmap to use
  for (i in 1:1000){
    print(paste("Now starting SIMMAP: ", i))
    tree <- func_simmaps[[i]]
    
    # apply function to each model
    results <- lapply(models, OUwie.model, phy=tree, data=OUwie_data) 
    q = 1
    
    for (q in 1:2) {
      each_model <- results[[q]] 
      
      # write results into data frame
      best_model_data_2$trait[r] <- trait_vector[j]
      best_model_data_2$simmap_count[r] <- i
      best_model_data_2$model[r] <- each_model$model
      best_model_data_2$loglik[r] <- each_model$loglik
      best_model_data_2$aicc[r] <- each_model$AICc
      best_model_data_2$eigval[r] <- list(each_model$eigval)
      best_model_data_2$sigma.sq[r] <- each_model$solution[2,1]
      best_model_data_2$alpha[r] <- each_model$solution[1,1] # keep?
      best_model_data_2$theta[r] <- each_model$theta[1,1]
      best_model_data_2$theta_se[r] <- each_model$theta[1,2]
      best_model_data_2$saddle[r] <- any(each_model$eigval < 0)
      
      q <- q + 1
      r <- r + 1
    }
    i <- i + 1
  }
  j <- j + 1
}

BM1_OU1<-best_model_data_2
save(BM1_OU1, file="BM1_OU1.RData")


#################################################################
#   OUwie - OUM BMS
#################################

###################  Set up data frame for results #############################

best_model_data <- tibble(.rows = 10000) # 5 traits * 2 models * 1000 simmaps
best_model_data$trait <- rep(NA,10000)
best_model_data$simmap_count <- rep(NA,10000)
best_model_data$model <- rep(NA,10000)
best_model_data$loglik <- rep(NA,10000)
best_model_data$aicc <- rep(NA,10000)
best_model_data$eigval <- rep(NA, 10000)
best_model_data$sigma.sq_ins <- rep(NA, 10000)
best_model_data$sigma.sq_inv <- rep(NA, 10000)
best_model_data$sigma.sq_carn <- rep(NA, 10000)
best_model_data$sigma.sq_hd <- rep(NA, 10000)
best_model_data$sigma.sq_larv <- rep(NA, 10000)
best_model_data$sigma.sq_omn <- rep(NA, 10000)
best_model_data$alpha <- rep(NA,10000)
best_model_data$theta_ins <- rep(NA,10000)
best_model_data$theta_inv <- rep(NA,10000)
best_model_data$theta_carn <- rep(NA,10000)
best_model_data$theta_hd <- rep(NA,10000)
best_model_data$theta_larv <- rep(NA,10000)
best_model_data$theta_omn <- rep(NA,10000)
best_model_data$theta_ins_se <- rep(NA,10000)
best_model_data$theta_inv_se <- rep(NA,10000)
best_model_data$theta_carn_se <- rep(NA,10000)
best_model_data$theta_hd_se <- rep(NA,10000)
best_model_data$theta_larv_se<- rep(NA,10000)
best_model_data$theta_omn_se <- rep(NA,10000)
best_model_data$saddle <- rep(NA,10000)


#########################  Run OUwie  ##########################################

library(OUwie)

# define function
OUwie.model <- function(model, phy, data) {
  print(paste("Now starting model: ", model))
  return(OUwie(phy, data, model, simmap.tree=T, diagn=T)) 
}

models <- c("OUM", "BMS")

trait_vector <- c("IL_log","IL_resid","ABV_log","ABV_resid", "Size")

j = 1 #trait
i = 1 #simmap
r = 1 
q = 1 #model

# subset data frame for each trait
for (j in 1:5){
  print(paste("Now starting trait: ", trait_vector[j]))
  OUwie_data <- data.frame(dat[ ,c(1,2,j+2)])  
  i = 1
  
  # extract simmap to use
  for (i in 1:1000){
    print(paste("Now starting SIMMAP: ", i))
    tree <- func_simmaps[[i]]
    
    # apply function to each model
    results <- lapply(models, OUwie.model, phy=tree, data=OUwie_data) 
    q = 1
    
    for (q in 1:2) {
      each_model <- results[[q]] 
      
      # write results into data frame
      best_model_data$trait[r] <- trait_vector[j]
      best_model_data$simmap_count[r] <- i
      best_model_data$model[r] <- each_model$model
      best_model_data$loglik[r] <- each_model$loglik
      best_model_data$aicc[r] <- each_model$AICc
      best_model_data$eigval[r] <- list(each_model$eigval)
      best_model_data$sigma.sq_ins[r] <- each_model$solution[2,1]
      best_model_data$sigma.sq_inv[r] <- each_model$solution[2,2]
      best_model_data$sigma.sq_carn[r] <- each_model$solution[2,3]
      best_model_data$sigma.sq_hd[r] <- each_model$solution[2,4]
      best_model_data$sigma.sq_larv[r] <- each_model$solution[2,5]
      best_model_data$sigma.sq_omn[r] <- each_model$solution[2,6]
      best_model_data$alpha[r] <- each_model$solution[1,1]
      best_model_data$theta_ins[r] <- each_model$theta[1,1]
      best_model_data$theta_inv[r] <- each_model$theta[2,1]
      best_model_data$theta_carn[r] <- each_model$theta[3,1]
      best_model_data$theta_hd[r] <- each_model$theta[4,1]
      best_model_data$theta_larv[r] <- each_model$theta[5,1]
      best_model_data$theta_omn[r] <- each_model$theta[6,1]
      best_model_data$theta_ins_se[r] <- each_model$theta[1,2]
      best_model_data$theta_inv_se[r] <- each_model$theta[2,2]
      best_model_data$theta_carn_se[r] <- each_model$theta[3,2]
      best_model_data$theta_hd_se[r] <- each_model$theta[4,2]
      best_model_data$theta_larv_se[r] <- each_model$theta[5,2]
      best_model_data$theta_omn_se[r] <- each_model$theta[6,2]
      best_model_data$saddle[r] <- any(each_model$eigval < 0)
      
      q<- q + 1
      r <- r + 1
    }
    i <- i + 1
  }
  j <- j + 1
}

OUM_BMS<-best_model_data
save(OUM_BMS, file="OUM_BMS.RData")




load("OUM_BMS.RData")
load("OU1_BM1.RData")

load("BM1_OU1_body size.RData")
load("OUM_BMS_function_body size.RData")



trait_vector <- c("IL_log","IL_resid","ABV_log","ABV_resid", "Size")


aicc_results <- setNames(data.frame(matrix(ncol = 5, nrow = 5)), c("trait","BM1", "BMS", "OU1", "OUM")) 


for(i in 1:length(trait_vector)){
  aicc_results$trait[i]<-trait_vector[i]
  OU1<-filter(BM1_OU1, trait == trait_vector[i], model=="OU1")
  aicc_results$OU1[i]<- mean(OU1$aicc)
  BM1<-filter(BM1_OU1, trait == trait_vector[i], model=="BM1")
  aicc_results$BM1[i]<-mean(BM1$aicc)
  OUM<-filter(OUM_BMS, trait == trait_vector[i], model=="OUM")
  aicc_results$OUM[i]<-mean(OUM$aicc)
  BMS<-filter(OUM_BMS, trait == trait_vector[i], model=="BMS")
  aicc_results$BMS[i]<-mean(BMS$aicc)
  
}

write.csv(aicc_results, file="AIcc_results.csv")

aicc_results_size <- setNames(data.frame(matrix(ncol = 5, nrow = 1)), c("trait","BM1", "BMS", "OU1", "OUM")) 

for(i in 1:1){
  aicc_results_size$trait[i]<- "Max Body Size"
  OU1<-filter(BM1_OU1_size, model=="OU1")
  aicc_results_size$OU1[i]<- mean(OU1$aicc)
  BM1<-filter(BM1_OU1_size, model=="BM1")
  aicc_results_size$BM1[i]<-mean(BM1$aicc)
  OUM<-filter(OUM_BMS_size, model=="OUM")
  aicc_results_size$OUM[i]<-mean(OUM$aicc)
  BMS<-filter(OUM_BMS_size,  model=="BMS")
  aicc_results_size$BMS[i]<-mean(BMS$aicc)
  
}


OUM.il<-filter(OUM_BMS, trait == "IL_log", model=="OUM")
theta_results_gape <- setNames(data.frame(matrix(ncol = 2, nrow = 6000)), c("Diet", "Theta")) 
theta_results_gape$Function[1:1000] <-"AquaticInsectivore"
theta_results_gape$Theta[1:1000] <-OUM.il$theta_ins
theta_results_gape$Function[1001:2000] <-"AquaticInvertivore"
theta_results_gape$Theta[1001:2000] <-OUM.il$theta_inv
theta_results_gape$Function[2001:3000] <-"Carnivore"
theta_results_gape$Theta[2001:3000] <-OUM.il$theta_carn
theta_results_gape$Function[3001:4000] <-"HD"
theta_results_gape$Theta[3001:4000] <-OUM.il$theta_hd
theta_results_gape$Function[4001:5000] <-"InsectLarvaphage"
theta_results_gape$Theta[4001:5000] <-OUM.il$theta_larv
theta_results_gape$Function[5001:6000] <-"Omnivore"
theta_results_gape$Theta[5001:6000] <-OUM.il$theta_omn


plot1<-ggplot(theta_results_gape, aes(x=Theta, fill=Function)) + 
  geom_density() +
  labs(title="Log Intestine Length") +
  theme(plot.title = element_text(hjust = 0.5))

OUM.il.resid<-filter(OUM_BMS, trait == "IL_resid", model=="OUM")
theta_results_gape <- setNames(data.frame(matrix(ncol = 2, nrow = 6000)), c("Diet", "Theta")) 
theta_results_gape$Function[1:1000] <-"AquaticInsectivore"
theta_results_gape$Theta[1:1000] <-OUM.il.resid$theta_ins
theta_results_gape$Function[1001:2000] <-"AquaticInvertivore"
theta_results_gape$Theta[1001:2000] <-OUM.il.resid$theta_inv
theta_results_gape$Function[2001:3000] <-"Carnivore"
theta_results_gape$Theta[2001:3000] <-OUM.il.resid$theta_carn
theta_results_gape$Function[3001:4000] <-"HD"
theta_results_gape$Theta[3001:4000] <-OUM.il.resid$theta_hd
theta_results_gape$Function[4001:5000] <-"InsectLarvaphage"
theta_results_gape$Theta[4001:5000] <-OUM.il.resid$theta_larv
theta_results_gape$Function[5001:6000] <-"Omnivore"
theta_results_gape$Theta[5001:6000] <-OUM.il.resid$theta_omn


plot2<-ggplot(theta_results_gape, aes(x=Theta, fill=Function)) + 
  geom_density() +
  labs(title="Size Corrected Intestine Length") +
  theme(plot.title = element_text(hjust = 0.5))


OUM.abv<-filter(OUM_BMS, trait == "ABV_log", model=="OUM")
theta_results_gape <- setNames(data.frame(matrix(ncol = 2, nrow = 6000)), c("Diet", "Theta")) 
theta_results_gape$Function[1:1000] <-"AquaticInsectivore"
theta_results_gape$Theta[1:1000] <-OUM.abv$theta_ins
theta_results_gape$Function[1001:2000] <-"AquaticInvertivore"
theta_results_gape$Theta[1001:2000] <-OUM.abv$theta_inv
theta_results_gape$Function[2001:3000] <-"Carnivore"
theta_results_gape$Theta[2001:3000] <-OUM.abv$theta_carn
theta_results_gape$Function[3001:4000] <-"HD"
theta_results_gape$Theta[3001:4000] <-OUM.abv$theta_hd
theta_results_gape$Function[4001:5000] <-"InsectLarvaphage"
theta_results_gape$Theta[4001:5000] <-OUM.abv$theta_larv
theta_results_gape$Function[5001:6000] <-"Omnivore"
theta_results_gape$Theta[5001:6000] <-OUM.abv$theta_omn


plot3<-ggplot(theta_results_gape, aes(x=Theta, fill=Function)) + 
  geom_density() +
  labs(title="Abdominal Cavity Volume") + xlim(2,5) +
  theme(plot.title = element_text(hjust = 0.5))


require(gridExtra)

pdf(file="PLot of theta values_final.pdf", height=8, width=8)
grid.arrange(plot1, plot2, plot3, ncol=1, nrow=3)
dev.off()


#########################################################
#         OUWIE sim
#########################################################

dd<-ReorderData(tree, dat, taxa.names="Species")

phy <- pbtree(n = 20, nsim = 2) 
disc_trait <- setNames(sample(letters[1:2], 20, replace = TRUE), phy[[1]]$tip.label)
cont_traits <- as_tibble(iris[1:20, 1:2]) %>%
  mutate(species = phy[[1]]$tip.label)

models <- c("BM1", "BMS", "OU1", "OUM") 

results <- ouwie_tidy(phy, disc_trait, cont_traits, models, nsim = 1)

i=1
OUwie_sim <- setNames(data.frame(matrix(ncol = 3, nrow = 123)), c("trait","BM1", "BMS", "OU1", "OUM")) 
for(i in 1:length(func_simmaps)){
OUwie
sim<-OUwie.sim(phy=func_simmaps[[i]],  simmap.tree=TRUE, root.age=NULL, scaleHeight=FALSE,
          alpha=c(0.8792516, 0.8792516, 0.8792516, 0.8792516, 0.8792516, 0.8792516), sigma.sq=c(0.0652916,0.0652916,0.0652916,0.0652916,0.0652916,0.0652916), theta0 = 0, 
          theta=c(1.570477, 1.695760, 1.882693, 2.336379, 1.674587, 1.663936), mserr="none", shift.point=0.5, 
          fitted.object=NULL, get.all=FALSE)

dd.sim<-ReorderData(tree, sim, taxa.names="Genus_species")


