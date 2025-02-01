require(fishtree)
require(phytools)
require(geiger)
require(ape)
require(geomorph)

tree<-fishtree_phylogeny()
dat3 <- readxl::read_xlsx("C:/Users/Test 1/OneDrive/Documents/AbCav Project/data_ouwie.xlsx")
rownames(dat3) <- dat3$Species
d<-name.check(tree, dat3) 

d
#tree <- read.tree("C:/Users/gabmc/Downloads/Schonhuth et al 2018_Leuciscidae/Schonhuth et al 2018_Leuciscidae/Leuciscid trees/RAxML_bestTree.LeuCombined_ML copy.tre")
phy.mcc <- drop.tip(tree,d$tree_not_data)
plot(phy.mcc, cex = 0.8)
write.nexus(phy.mcc, file = "C:/Users/Test 1/OneDrive/Documents/AbCav Project/droptiptree.nexus", translate = TRUE)

#testing testing123
#par(mfrow=c(1,1))
#plotTree(phy.mcc,ftype="i",lwd=1,fsize=0.6,type="fan",part=0.88, 
    #     tip.color="green")


#making data getting rid of species that weren't in the tree
#dat2 <- dat[rownames(dat) %in% phy.mcc$tip.label,] 
#rownames(dat2) <- dat2$Species
#size<-dat$`Average Length`
#names(size)<-rownames(dat)


#labelling size and maximum length variables + trophic categories
size<-dat3$`Standard Length`
names(size)<-rownames(dat3)
maxlgn<-dat3$`Max Length`
names(maxlgn)<-rownames(dat3)
diet<-dat3$Trophic
names(diet)<-rownames(dat3)
acv <- dat3$ABV
names(acv) <- rownames(dat3)
itl <- dat3$`Intestine Length`
names(itl) <- rownames(dat3)
acv.log <- dat3$ABV_log
names(acv.log) <- rownames(dat3)
itl.log <- dat3$IL_log
names(itl.log) <- rownames(dat3)
sc.itl <- dat3$IL_resid
names(sc.itl) <- rownames(dat3)
sc.acv <- dat3$ABV_resid
names(sc.acv) <- rownames(dat3)

#cv<-phyl.resid(phy.mcc, log(size), log(cav.vol))
#il<-phyl.resid(phy.mcc, log(size), log(int.lgn))
#cv$resid
#il$resid

#saving variables as separate files
#write.csv(diet, "C:/Users/Test 1/OneDrive/Documents/AbCav Project/trophicgroup.csv")
#write.csv(size, "C:/Users/Test 1/OneDrive/Documents/AbCav Project/measuredlength.csv")
#write.csv(maxlgn, "C:/Users/Test 1/OneDrive/Documents/AbCav Project/maxlength.csv")
#write.csv(int.lgn, "C:/Users/Test 1/OneDrive/Documents/AbCav Project/intestinelength.csv")
#write.csv(cav.vol, "C:/Users/Test 1/OneDrive/Documents/AbCav Project/abcavvolume.csv")
#write.csv(cv$resid, "C:/Users/Test 1/OneDrive/Documents/AbCav Project/cavvolresid.csv")
#write.csv(il$resid, "C:/Users/Test 1/OneDrive/Documents/AbCav Project/intlgnresid.csv")

#graph display (columns, rows)
par(mfrow=c(3,1))

#plot body size vs intestine length
plot(log(size), itl.log, 
     main = "Body Length vs Intestine Length",
     xlab = "Log-transformed Body Length (mm)",
     ylab = "Log-transformed Intestine Length (mm)",
     pch = 19,
     text(3.6, 2.65, labels="Rsq=0.647"),
     col = ifelse(dat2$`Trophic/Diet Guild` == "HD", "#33a02c", 
                         ifelse(dat2$`Trophic/Diet Guild` == "AquaticInvertivore", "#b2df8a", 
                                ifelse(dat2$`Trophic/Diet Guild` == "Omnivore", "grey", 
                                       ifelse(dat2$`Trophic/Diet Guild` == "AquaticInsectivore", "#a6cee3",
                                              ifelse(dat2$`Trophic/Diet Guild` == "InsectLarvaphage", "#1f78b4", "#363636"))))))

legend(x = "topleft", title="Legend", cex = 0.6,
  legend=c("Herbivore-Detritivore", "Aquatic Invertivore", "Aquatic Insectivore", "Insect Larvaphage", "Omnivore", "Carnivore"), 
  fill = c("#33a02c", "#b2df8a", "#a6cee3", "#1f78b4", "grey", "#363636"))
pglsmodel1<-procD.pgls(log(size) ~ itl.log, phy=phy.mcc, iter = 10000)
summary(pglsmodel1)
plot(pglsmodel1, type="regression", reg.type="RegScore", predictor = dat3$IL_log)
abline(pglsmodel1)
#abline(lm(itl.log~log(size)))
?abline

#log adjusted of body length vs abdominal cavity:
plot(log(size), acv.log, 
     main = "Body Length vs Abdominal Cavity Volume",
     xlab = "Log-transformed Body Length (mm)",
     ylab = "Log-transformed Abdominal Cavity Volume (mm3)",
     pch = 19,
     text(3.6, 3.75, labels="Rsq=0.878"),
     col = ifelse(dat2$`Trophic/Diet Guild` == "HD", "#33a02c", 
                  ifelse(dat2$`Trophic/Diet Guild` == "AquaticInvertivore", "#b2df8a", 
                         ifelse(dat2$`Trophic/Diet Guild` == "Omnivore", "grey", 
                                ifelse(dat2$`Trophic/Diet Guild` == "AquaticInsectivore", "#a6cee3",
                                       ifelse(dat2$`Trophic/Diet Guild` == "InsectLarvaphage", "#1f78b4", "#363636"))))))

legend(x = "topleft", title="Legend", cex = 0.6,
       legend=c("Herbivore-Detritivore", "Aquatic Invertivore", "Aquatic Insectivore", "Insect Larvaphage", "Omnivore", "Carnivore"), 
       fill = c("#33a02c", "#b2df8a", "#a6cee3", "#1f78b4", "grey", "#363636"))
pglsmodel2<-procD.pgls(log(size) ~ acv.log, phy=phy.mcc, iter = 10000)
summary(pglsmodel2)
abline(pglsmodel2)
#abline(lm(acv.log~log(size)))

#log adjusted of intestine length vs abdominal cavity:
plot(itl.log, acv.log, 
     main = "Intestine Length vs Abdominal Cavity Volume",
     xlab = "Log-transformed Intestine Length (mm)",
     ylab = "Log-transformed Abdominal Cavity Volume (mm3)",
     pch = 19,
     text(1.35, 3.75, labels="Rsq=0.489"),
     col = ifelse(dat2$`Trophic/Diet Guild` == "HD", "#33a02c", 
                  ifelse(dat2$`Trophic/Diet Guild` == "AquaticInvertivore", "#b2df8a", 
                         ifelse(dat2$`Trophic/Diet Guild` == "Omnivore", "grey", 
                                ifelse(dat2$`Trophic/Diet Guild` == "AquaticInsectivore", "#a6cee3",
                                       ifelse(dat2$`Trophic/Diet Guild` == "InsectLarvaphage", "#1f78b4", "#363636"))))))

legend(x = "topleft", title="Legend", cex = 0.6,
       legend=c("Herbivore-Detritivore", "Aquatic Invertivore", "Aquatic Insectivore", "Insect Larvaphage", "Omnivore", "Carnivore"), 
       fill = c("#33a02c", "#b2df8a", "#a6cee3", "#1f78b4", "grey", "#363636"))
pglsmodel3<-procD.pgls(itl.log ~ acv.log, phy=phy.mcc, iter = 10000)
summary(pglsmodel3)
abline(pglsmodel3)
#abline(lm(acv.log~itl.log))

#graph display (columns, rows)
par(mfrow=c(1,1))

#legend lol
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
#legend("topleft", title="Legend", legend =c("Aquatic Insectivore", "Aquatic Invertivore", "Carnivore", "Herbivore/Detritivore", "Insect Larvaphage", "Omnivore"), 
       #fill = c("#E69F00", "#0072B2", "grey", "#009E73", "#CC79A7", "#56B4E9"))
legend('top', fill=c("#33a02c", "#b2df8a", "#a6cee3", "#1f78b4", "grey", "#363636"), legend=c("Herbivore-Detritivore", "Aq. Invertivore", "Aq. Insectivore", "Insect Larvaphage", "Omnivore", "Carnivore"),
     inset=c(0, -.15), xpd=TRUE, ncol=3)
#legend(x = "top", horiz=TRUE, inset=c(0, 0), xpd=TRUE,
      # legend=c("Herbivore-Detritivore", "Aq. Invertivore", "Aq. Insectivore", "Insect Larvaphage", "Omnivore", "Carnivore"), 
      # fill = c("#33a02c", "#b2df8a", "#a6cee3", "#1f78b4", "grey", "#363636"))


##BOXPLOTS

par(mfrow=c(1,3))
#size corrected int lgn by trophic group
boxplot(sc.itl ~ dat3$Trophic,
        xaxt = "n",
        data = dat3,
        main = "Intestine Length",
        xlab = "Trophic Group",
        ylab = "Size-Corrected Intestine Length",
        col = c("#a6cee3", "#b2df8a", "#363636", "#33a02c", "#1f78b4", "grey"),
        names = c("Aq. Insectivore","Aq. Invertivore","Carnivore", "HD", "Insect Larvaphage", "Omnivore"))
legend(x = "bottom", cex = 1, horiz=TRUE, inset=c(0, -.15), xpd=TRUE,
       legend=c("Herbivore-Detritivore", "Aquatic Invertivore", "Aquatic Insectivore", "Insect Larvaphage", "Omnivore", "Carnivore"), 
       fill = c("#33a02c", "#b2df8a", "#a6cee3", "#1f78b4", "grey", "#363636"))
eco<-dat3$Trophic
names(eco)<-rownames(dat3)
trophic.anova<-procD.pgls(sc.itl~eco, phy=phy.mcc, iter = 10000)
summary(trophic.anova)
gp <-  interaction (eco)
PW <- pairwise(trophic.anova, groups = gp, covariate = NULL)
summary(PW)

#size-corrected ab cav by trophic group
boxplot(sc.acv ~ dat3$Trophic,
        xaxt = "n",
        data = dat3,
        main = "Abdominal Cavity Volume",
        xlab = "Trophic Group",
        ylab = "Size-Corrected Abdominal Cavity Volume",
        col = c("#a6cee3", "#b2df8a", "#363636", "#33a02c", "#1f78b4", "grey"),
        names = c("Aq. Insectivore","Aq. Invertivore","Carnivore", "HD", "Insect Larvaphage", "Omnivore"))
legend(x = "topleft", title="Legend", cex = 0.6,
       legend=c("Herbivore-Detritivore", "Aquatic Invertivore", "Aquatic Insectivore", "Insect Larvaphage", "Omnivore", "Carnivore"), 
       fill = c("#33a02c", "#b2df8a", "#a6cee3", "#1f78b4", "grey", "#363636"))
eco<-dat3$Trophic
names(eco)<-rownames(dat3)
trophic.anova<-procD.pgls(sc.acv~eco, phy=phy.mcc, iter = 10000)
summary(trophic.anova)
gp <-  interaction (eco)
PW <- pairwise(trophic.anova, groups = gp, covariate = NULL)
summary(PW)

#max body length by trophic group
boxplot(log(maxlgn) ~ dat3$Trophic,
        xaxt = "n",
        data = dat3,
        main = "Maximum Body Length",
        xlab = "Trophic Group",
        ylab = "Log-adjusted Maximum Body Length (mm)",
        col = c("#a6cee3", "#b2df8a", "#363636", "#33a02c", "#1f78b4", "grey"),
        names = c("Aq. Insectivore","Aq. Invertivore","Carnivore", "HD", "Insect Larvaphage", "Omnivore"))
legend(x = "topleft", title="Legend", cex = 0.6,
       legend=c("Herbivore-Detritivore", "Aquatic Invertivore", "Aquatic Insectivore", "Insect Larvaphage", "Omnivore", "Carnivore"), 
       fill = c("#33a02c", "#b2df8a", "#a6cee3", "#1f78b4", "grey", "#363636"))
eco<-dat3$Trophic
names(eco)<-rownames(dat3)
trophic.anova<-procD.pgls(log(maxlgn)~eco, phy=phy.mcc, iter = 10000)
summary(trophic.anova)
gp <-  interaction (eco)
PW <- pairwise(trophic.anova, groups = gp, covariate = NULL)
summary(PW)

#par(mfrow=c(1,1))
#measured body length by trophic group
boxplot(log(size) ~ dat3$Trophic,
        xaxt = "n",
        data = dat3,
        main = "Body Length by Trophic Group",
        xlab = "Trophic Group",
        ylab = "Body Length Log",
        col = c("#a6cee3", "#b2df8a", "black", "#33a02c", "#1f78b4", "grey"),
        names = c("Aq. Insectivore","Aq. Invertivore","Carnivore", "HD", "Insect Larvaphage", "Omnivore"))
legend(x = "topleft", title="Legend", cex = 0.6,
       legend=c("Herbivore-Detritivore", "Aquatic Invertivore", "Aquatic Insectivore", "Insect Larvaphage", "Omnivore", "Carnivore"), 
       fill = c("#33a02c", "#b2df8a", "#a6cee3", "#1f78b4", "grey", "black"))
eco<-dat3$Trophic
names(eco)<-rownames(dat3)
trophic.anova<-procD.pgls(log(size)~eco, phy=phy.mcc, iter = 10000)
summary(trophic.anova)
gp <-  interaction (eco)
PW <- pairwise(trophic.anova, groups = gp, covariate = NULL)
summary(PW)


#playing around with color
plot(log(size), itl.log, 
     main = "Body Length vs Intestine Length",
     xlab = "Log-transformed Body Length (mm)",
     ylab = "Log-transformed Intestine Length (mm)",
     pch = 19,
     col = ifelse(dat2$`Trophic/Diet Guild` == "HD", "#33a02c", 
                         ifelse(dat2$`Trophic/Diet Guild` == "AquaticInvertivore", "#b2df8a", 
                                ifelse(dat2$`Trophic/Diet Guild` == "Omnivore", "grey", 
                                       ifelse(dat2$`Trophic/Diet Guild` == "AquaticInsectivore", "#a6cee3",
                                              ifelse(dat2$`Trophic/Diet Guild` == "InsectLarvaphage", "#1f78b4", "black"))))))

legend(x = "topleft", title="Legend", cex = 0.75,
       legend=c("Aquatic Insectivore", "Aquatic Invertivore", "Carnivore", "Herbivore/Detritivore", "Insect Larvaphage", "Omnivore"), 
       fill = c("#a6cee3", "#b2df8a", "black", "#33a02c", "#1f78b4", "grey"))

plot(log(size), itl.log, 
     main = "Body Length vs Intestine Length",
     xlab = "Log-transformed Body Length (mm)",
     ylab = "Log-transformed Intestine Length (mm)",
     pch = 19,
     col = ifelse(dat2$`Trophic/Diet Guild` == "HD", "#1b7837", 
                  ifelse(dat2$`Trophic/Diet Guild` == "AquaticInvertivore", "#7fbf7b", 
                         ifelse(dat2$`Trophic/Diet Guild` == "Omnivore", "#af8dc3", 
                                ifelse(dat2$`Trophic/Diet Guild` == "AquaticInsectivore", "#d9f0d3",
                                       ifelse(dat2$`Trophic/Diet Guild` == "InsectLarvaphage", "#e7d4e8", "#762a83"))))))

legend(x = "topleft", title="Legend", cex = 0.75,
       legend=c("Aquatic Insectivore", "Aquatic Invertivore", "Carnivore", "Herbivore/Detritivore", "Insect Larvaphage", "Omnivore"), 
       fill = c("#d9f0d3", "#7fbf7b", "#762a83", "#1b7837", "#e7d4e8", "#af8dc3"))


######################supplemental regression plots
par(mfrow=c(3,1))

#plot body size vs intestine length
plot(log(size), sc.itl, 
     main = "Body Length vs Size-Corrected Intestine Length",
     xlab = "Log-transformed Body Length (mm)",
     ylab = "Size-Corrected Intestine Length (mm)",
     pch = 19,
     text(3.6, 1.2, labels="Rsq=0"),
     col = ifelse(dat2$`Trophic/Diet Guild` == "HD", "#33a02c", 
                  ifelse(dat2$`Trophic/Diet Guild` == "AquaticInvertivore", "#b2df8a", 
                         ifelse(dat2$`Trophic/Diet Guild` == "Omnivore", "grey", 
                                ifelse(dat2$`Trophic/Diet Guild` == "AquaticInsectivore", "#a6cee3",
                                       ifelse(dat2$`Trophic/Diet Guild` == "InsectLarvaphage", "#1f78b4", "#363636"))))))

legend(x = "topleft", title="Legend", cex = 0.6,
       legend=c("Herbivore-Detritivore", "Aquatic Invertivore", "Aquatic Insectivore", "Insect Larvaphage", "Omnivore", "Carnivore"), 
       fill = c("#33a02c", "#b2df8a", "#a6cee3", "#1f78b4", "grey", "#363636"))
pglsmodel4<-procD.pgls(log(size) ~ sc.itl, phy=phy.mcc, iter = 10000)
summary(pglsmodel4)
abline(pglsmodel)
#plot(pglsmodel1, type="regression", reg.type="RegScore", predictor = dat3$IL_log)
#abline(pglsmodel1)
#abline(lm(itl.log~log(size)))
#?abline

#log adjusted of body length vs abdominal cavity:
plot(log(size), sc.acv, 
     main = "Body Length vs Size-Corrected Abdominal Cavity Volume",
     xlab = "Log-transformed Body Length (mm)",
     ylab = "Size-Corrected Abdominal Cavity Volume (mm3)",
     pch = 19,
     text(3.6, 0.7, labels= "Rsq=0"),
     col = ifelse(dat2$`Trophic/Diet Guild` == "HD", "#33a02c", 
                  ifelse(dat2$`Trophic/Diet Guild` == "AquaticInvertivore", "#b2df8a", 
                         ifelse(dat2$`Trophic/Diet Guild` == "Omnivore", "grey", 
                                ifelse(dat2$`Trophic/Diet Guild` == "AquaticInsectivore", "#a6cee3",
                                       ifelse(dat2$`Trophic/Diet Guild` == "InsectLarvaphage", "#1f78b4", "#363636"))))))

legend(x = "topleft", title="Legend", cex = 0.6,
       legend=c("Herbivore-Detritivore", "Aquatic Invertivore", "Aquatic Insectivore", "Insect Larvaphage", "Omnivore", "Carnivore"), 
       fill = c("#33a02c", "#b2df8a", "#a6cee3", "#1f78b4", "grey", "#363636"))
pglsmodel5<-procD.pgls(log(size) ~ sc.acv, phy=phy.mcc, iter = 10000)
summary(pglsmodel5)
abline(pglsmodel5)
#abline(lm(acv.log~log(size)))

#log adjusted of intestine length vs abdominal cavity:
plot(sc.itl, sc.acv, 
     main = "Size-Corrected Intestine Length vs Size-Corrected Abdominal Cavity Volume",
     xlab = "Size-Corrected Intestine Length",
     ylab = "Size-Corrected Abdominal Cavity Volume",
     pch = 19,
     text(-0.85, 0.7, labels="R=0.06"),
     col = ifelse(dat2$`Trophic/Diet Guild` == "HD", "#33a02c", 
                  ifelse(dat2$`Trophic/Diet Guild` == "AquaticInvertivore", "#b2df8a", 
                         ifelse(dat2$`Trophic/Diet Guild` == "Omnivore", "grey", 
                                ifelse(dat2$`Trophic/Diet Guild` == "AquaticInsectivore", "#a6cee3",
                                       ifelse(dat2$`Trophic/Diet Guild` == "InsectLarvaphage", "#1f78b4", "#363636"))))))

legend(x = "topleft", title="Legend", cex = 0.6,
       legend=c("Herbivore-Detritivore", "Aquatic Invertivore", "Aquatic Insectivore", "Insect Larvaphage", "Omnivore", "Carnivore"), 
       fill = c("#33a02c", "#b2df8a", "#a6cee3", "#1f78b4", "grey", "#363636"))
pglsmodel6<-procD.pgls(sc.itl ~ sc.acv, phy=phy.mcc, iter = 10000)
summary(pglsmodel6)
abline(pglsmodel6)
#testing

reg_figure <- ggarrange(reg1, reg2, reg3,
                        labels = c("a", "b", "c"),
                        ncol = 1, nrow = 3)
reg_figure

box_figure <- ggarrange(box1, box2, box3,
                        labels = c("a", "b", "c"),
                        ncol = 3, nrow = 1)
box_figure

sup_figure

#testing git/github