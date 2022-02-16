###########################################
######### PCA AND DAPC ANALYSES ###########
########## By: Chan Kin Onn ############### 
###########################################

library(ggplot2)
setwd("~/Desktop/Pubisulcus_complex/morphology/datasets/mar2/3_region/pca/")
pcadata <- read.csv("3region_adj_log.csv") #read table that is in csv format
pcadata2 <- log10(pcadata[,2:ncol(pcadata)]) #log-transform pcadata
species <- pcadata[,1] #isolate colum with OTU's
pcadata3 <- cbind(species, pcadata2) #create new pcadata frame with OTU's and log-transformed pcadata
pca <- prcomp(pcadata3[,2:ncol(pcadata3)], scale=TRUE) #perform PCA with scaling. 
##Scaling uses the std dev as a scaling factor. After scaling, all pcadata have a std dev of one, therefore pcadata is
##analyzed on the basis of correlations not covariances, as is the case with centering.
##variables in different units or same units but have large variances must be scaled.
scores <- data.frame(species, pca$x[,1:10]) #create pcadata frame using species and PCA scores
str(scores)

#Set colors
cbPalette <- c("#999999", "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #Colorblind friendly
myCol <- c("#009E73", "#56B4E9", "#E69F00") #Custom color pallette

a <- ggplot(scores, aes(x=PC1, y=PC2, group = pcadata3$species)) +
  theme_bw() + theme(legend.key = element_blank()) +
  scale_color_manual(values = myCol) + stat_ellipse() +
  geom_point((aes(color=pcadata$clade)), size=6) + xlab("PC1 37%") + ylab("PC2 11%")+
  theme(legend.title=element_blank(), axis.text=element_text(size=16), axis.title=element_text(size=16),
        legend.text=element_text(size=24))+guides(colour = guide_legend(override.aes = list(size=6)))
a

b <- ggplot(scores, aes(x=PC1, y=PC3, group = pcadata3$species)) +
  theme_bw() + theme(legend.key = element_blank()) +
  scale_color_manual(values = myCol)  + stat_ellipse() +
  geom_point((aes(color=pcadata$clade)), size=6) + xlab("PC 11%") + ylab("PC3 9%")+
  theme(legend.title=element_blank(), axis.text=element_text(size=16), axis.title=element_text(size=16),
        legend.text=element_text(size=24))+guides(colour = guide_legend(override.aes = list(size=6)))
b

load <- pca$rotation
temp <- summary(pca) 
sum <- temp$importance #Extract important variables
eigen <- pca$sdev^2 #Extract eigenvalues
final_table <- rbind(sum, eigen, load) #Combine eigenvalues, loadings, and summary stats into single table
write.csv(sum, "4region_eigen.csv") # Enter your .csv file name.

ggsave(plot=b,height=10,width=10,dpi=300, filename="3region_PC2_3.pdf", useDingbats=FALSE)

## Plot 95% confidence ellipses
x <- pca$x[,1]
y <- pca$x[,2]
group <- pcadata[,1]
df <- data.frame(x=x, y=y, group=factor(group))
df_ell <- data.frame() 
for(g in levels(df$group)){df_ell <- rbind(df_ell, cbind(as.data.frame(with(df[df$group==g,], 
          ellipse(cor(x, y),scale=c(sd(x),sd(y)),centre=c(mean(x),mean(y))))),group=g))}

library(MASS) ##only if you want plot confidence ellipses##
library(ellipse) ##only if you want plot confidence ellipses##
library(devtools) ##only if you want plot confidence ellipses##
library(digest) ##only if you want plot confidence ellipses##
library(proto)##only if you want plot confidence ellipses##
library(gridExtra)

#Plot PCA scores/separate script
p <- ggplot(scores, aes(x=PC1, y=PC2, group=species)) + theme_bw() + theme(legend.key = element_blank()) + scale_color_manual(values=cbPalette) + scale_shape_manual(values=c(19)) + 
  geom_point(aes(color=species), size=6) + 
  theme(legend.title=element_blank(), axis.text=element_text(size=18), axis.title=element_text(size=20), 
        legend.text=element_text(size=18))+guides(colour = guide_legend(override.aes = list(size=8)))

p
#Write PCA results to table
loadings <- pca$rotation #Get PCA loadings
temp <- summary(pca) #Get PCA summary statistics. Summary() is not a pcadata frame or matrix. Cannot coerce. See next line to extract what we need
sum <- temp$importance #Extract important variables from summary
eigen <- pca$sdev^2 #Get eigenvalues
final.table <- rbind(sum, eigen, loadings) #Combine eigenvalues, loadings, and summary stats into single table
write.csv(final.table, "Borneo_limitedSZcorrect_eigen.csv") #Export table

##Export PDF for illustrator. Only for ggplot2
ggsave(plot=p,height=10,width=10,dpi=300, filename="Cyrtodactylus_Borneo.pdf", useDingbats=FALSE)

#Export JPG
jpeg("CyrtRatio_PCA.jpg", width=2200, height=2000, res=300)
ggplot(scores, aes(x=PC1, y=PC2, group=species)) + theme_bw() + theme(legend.key = element_blank()) + scale_shape_manual(values=c(1,3,5,6,7,8,9,10,4,2,19,15,16,17)) + 
   geom_point(aes(shape=species, color=species), size=6) + 
  theme(legend.title=element_blank(), axis.text=element_text(size=18), axis.title=element_text(size=20), 
        legend.text=element_text(size=15))+guides(colour = guide_legend(override.aes = list(size=8)))
dev.off()

###################################
######## DAPC analysis ############
###################################
library(adegenet)
dapc <- dapc(pcadata3[,2:ncol(pcadata)], grp=species)
scatter(dapc, col=cbPalette, scree.da=FALSE, scree.pca=FALSE, pch=1:13, clab=0, 
        cstar=0, leg=TRUE, cleg=1.4, cex=2.5, solid=0.85, posi.leg="bottomright") #cleg=size of legend

## Plot only the first discriminant function
scatter(dapc,1,1, col=myCol, scree.da=FALSE, clab=0, cstar=0, leg=TRUE, cleg=1.4, 
        cex=3.5, solid=0.85, posi.leg="topleft")

### Export PDF ###
pdf("dapc.female.nostar.pdf", height=7.5, width=9)
scatter(dapc, col=myCol, scree.da=FALSE, pch=15:17, clab=0, cstar=0, leg=TRUE, cleg=1.4, cex=2.5, solid=0.85, posi.leg="topleft")
dev.off()

