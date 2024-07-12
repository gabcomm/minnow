
require(phytools)
require(evobiR)
require(geiger)
require(OUwie)

tree<-read.nexus("tree.nexus")
dat<-read.csv("data_reduced.csv", row.names=1)

dd<-name.check(tree, dat, data.names=dat$Species)
tree.final<-drop.tip(tree,dd$tree_not_data)



#####################################################################
#               Simulations for IL resid
#####################################################################
dat.reorder<-ReorderData(tree.final, dat, taxa.names="row names")

sim.dat <- setNames(data.frame(matrix(ncol = 3, nrow = 107)), c("Species","Trophic", "Trait")) 

OUM.il.resid<-filter(OUM_BMS, trait == "IL_resid", model=="OUM")

sig<-mean(OUM.il.resid$sigma.sq_ins)
alp<-mean(OUM.il.resid$alpha)
ins<-mean(OUM.il.resid$theta_ins)
inv<-mean(OUM.il.resid$theta_inv)
carn<-mean(OUM.il.resid$theta_carn)
hd<-mean(OUM.il.resid$theta_hd)
larv<-mean(OUM.il.resid$theta_larv)
omn<-mean(OUM.il.resid$theta_omn)


sim_results_resid_IL <- setNames(data.frame(matrix(ncol = 6, nrow = 1000)), c("AquaticInsectivore","AquaticInvertivore", "Carnivore", "HD", "InsectLarvaphage", "Omnivore")) 

aic_results_resid_IL<-setNames(data.frame(matrix(ncol = 4, nrow = 1000)), c("BM1","BMS", "OU1", "OUM"))

for(i in 1:1000){
sim<-OUwie.sim(phy=func_simmaps[[i]],  simmap.tree=TRUE, root.age=NULL, scaleHeight=FALSE,
               alpha=c(alp, alp, alp, alp, alp, alp), sigma.sq=c(sig,sig,sig,sig,sig,sig), theta0 = 0, 
               theta=c(ins, inv, carn, hd, larv, omn), mserr="none", shift.point=0.5, 
               fitted.object=NULL, get.all=FALSE)

sim.dat$Species<-row.names(dat.reorder)
sim.dat$Trophic<-dat.reorder$Trophic

sim.dat$Trait<-sim$X +3

BM1<-OUwie(func_simmaps[[i]], sim.dat, model="BM1", simmap.tree=T, algorithm="three.point", diagn=T)
BMS<-OUwie(func_simmaps[[i]], sim.dat, model="BMS", simmap.tree=T, algorithm="three.point", diagn=T)
OU1<-OUwie(func_simmaps[[i]], sim.dat, model="OU1", simmap.tree=T, algorithm="three.point", diagn=T)
OUM<-OUwie(func_simmaps[[i]], sim.dat, model="OUM", simmap.tree=T, algorithm="three.point", diagn=T)

aic_results_resid_IL$BM1[i]<-BM1$AICc 
aic_results_resid_IL$BMS[i]<-BMS$AICc
aic_results_resid_IL$OU1[i]<-OU1$AICc
aic_results_resid_IL$OUM[i]<-OUM$AICc

sim_results_resid_IL$AquaticInsectivore[i]<-OUM$theta[1,1]
sim_results_resid_IL$AquaticInvertivore[i] <-OUM$theta[2,1]
sim_results_resid_IL$Carnivore[i] <-OUM$theta[3,1]
sim_results_resid_IL$HD[i] <-OUM$theta[4,1]
sim_results_resid_IL$InsectLarvaphage[i] <-OUM$theta[5,1]
sim_results_resid_IL$Omnivore[i] <-OUM$theta[6,1]

print(i)

}

write.csv(aic_results_resid_IL, file="AIC scores_il_resid.csv")
write.csv(sim_results_resid_IL, file="Theta values_il_resid.csv")
#####################################################################
#               Simulations for Log ABV
#####################################################################

OUM.abv<-filter(OUM_BMS, trait == "ABV_log", model=="OUM")


sim.dat <- setNames(data.frame(matrix(ncol = 3, nrow = 107)), c("Species","Trophic", "Trait")) 


sig<-mean(OUM.abv$sigma.sq_ins)
alp<-mean(OUM.abv$alpha)
ins<-mean(OUM.abv$theta_ins)
inv<-mean(OUM.abv$theta_inv)
carn<-mean(OUM.abv$theta_carn)
hd<-mean(OUM.abv$theta_hd)
larv<-mean(OUM.abv$theta_larv)
omn<-mean(OUM.abv$theta_omn)


sim_results_log_abv <- setNames(data.frame(matrix(ncol = 6, nrow = 1000)), c("AquaticInsectivore","AquaticInvertivore", "Carnivore", "HD", "InsectLarvaphage", "Omnivore")) 

aic_results_log_abv<-setNames(data.frame(matrix(ncol = 4, nrow = 1000)), c("BM1","BMS", "OU1", "OUM"))

for(i in 1:1000){
  sim<-OUwie.sim(phy=func_simmaps[[i]],  simmap.tree=TRUE, root.age=NULL, scaleHeight=FALSE,
                 alpha=c(alp, alp, alp, alp, alp, alp), sigma.sq=c(sig,sig,sig,sig,sig,sig), theta0 = 0, 
                 theta=c(ins, inv, carn, hd, larv, omn), mserr="none", shift.point=0.5, 
                 fitted.object=NULL, get.all=FALSE)
  
  sim.dat$Species<-row.names(dat.reorder)
  sim.dat$Trophic<-dat.reorder$Trophic
  
  sim.dat$Trait<-sim$X
  
  BM1<-OUwie(func_simmaps[[i]], sim.dat, model="BM1", simmap.tree=T, algorithm="three.point", diagn=T)
  BMS<-OUwie(func_simmaps[[i]], sim.dat, model="BMS", simmap.tree=T, algorithm="three.point", diagn=T)
  OU1<-OUwie(func_simmaps[[i]], sim.dat, model="OU1", simmap.tree=T, algorithm="three.point", diagn=T)
  OUM<-OUwie(func_simmaps[[i]], sim.dat, model="OUM", simmap.tree=T, algorithm="three.point", diagn=T)
  
  aic_results_log_abv$BM1[i]<-BM1$AICc 
  aic_results_log_abv$BMS[i]<-BMS$AICc
  aic_results_log_abv$OU1[i]<-OU1$AICc
  aic_results_log_abv$OUM[i]<-OUM$AICc
  
  sim_results_log_abv$AquaticInsectivore[i]<-OUM$theta[1,1]
  sim_results_log_abv$AquaticInvertivore[i] <-OUM$theta[2,1]
  sim_results_log_abv$Carnivore[i] <-OUM$theta[3,1]
  sim_results_log_abv$HD[i] <-OUM$theta[4,1]
  sim_results_log_abv$InsectLarvaphage[i] <-OUM$theta[5,1]
  sim_results_log_abv$Omnivore[i] <-OUM$theta[6,1]
  
  print(i)
  
}

write.csv(aic_results_log_abv, file="AIC scores_abv_log.csv")
write.csv(sim_results_log_abv, file="Theta values_abv_log.csv")
#####################################################################
#               Simulations for Resid ABV
#####################################################################

OUM.abv.resid<-filter(BM1_OU1, trait == "ABV_resid", model=="OU1")


sim.dat <- setNames(data.frame(matrix(ncol = 3, nrow = 107)), c("Species","Trophic", "Trait")) 


aic_results_resid_abv<-setNames(data.frame(matrix(ncol = 4, nrow = 1000)), c("BM1","BMS", "OU1", "OUM"))

for(i in 1:1000){
  sim<-OUwie.sim(phy=func_simmaps[[i]],  simmap.tree=TRUE, root.age=NULL, scaleHeight=FALSE,
                 alpha=c(0.577,0.577,0.577,0.577,0.577,0.577), sigma.sq=c(0.125,0.125,0.125,0.125,0.125,0.125), theta0 = 0, 
                 theta=c(-0.161,-0.161,-0.161,-0.161,-0.161,-0.161), mserr="none", shift.point=0.5, 
                 fitted.object=NULL, get.all=FALSE)
  
  sim.dat$Species<-row.names(dat.reorder)
  sim.dat$Trophic<-dat.reorder$Trophic
  
  sim.dat$Trait<-sim$X + 3
  
  BM1<-OUwie(func_simmaps[[i]], sim.dat, model="BM1", simmap.tree=T, algorithm="three.point", diagn=T)
  BMS<-OUwie(func_simmaps[[i]], sim.dat, model="BMS", simmap.tree=T, algorithm="three.point", diagn=T)
  OU1<-OUwie(func_simmaps[[i]], sim.dat, model="OU1", simmap.tree=T, algorithm="three.point", diagn=T)
  OUM<-OUwie(func_simmaps[[i]], sim.dat, model="OUM", simmap.tree=T, algorithm="three.point", diagn=T)
  
  aic_results_resid_abv$BM1[i]<-BM1$AICc 
  aic_results_resid_abv$BMS[i]<-BMS$AICc
  aic_results_resid_abv$OU1[i]<-OU1$AICc
  aic_results_resid_abv$OUM[i]<-OUM$AICc
  
  
  print(i)
  
}

write.csv(aic_results_resid_abv, file="AIC scores_abv_residual.csv")