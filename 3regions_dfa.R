library(MASS)
library(as.color)
library(klaR)
library(ape)
library(ggplot2)

setwd("~/Desktop/Pubisulcus_complex/Morphology/datasets/mar2/3_region/dfa/")
morpho <- read.csv(file= "3region_adj_log.csv", header=TRUE, fill=T, na.strings="")

cbPalette3 <- c("#999999", "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")

popcols <- morpho$clade
levels(popcols) <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#999999")
cols <- as.vector(popcols)

#try doing PCA on residuals in order to account for SVL
#try doing ANOVA on LD or ANOVA loadings.

##for quick plot 
# DFA by pop
# no training
head(morpho)
#data <- morpho[-c(1,3)]
data <- morpho[,2:ncol(morpho)] # REMOVE NA variables
data <- morpho[-c(1,3,12:14)] # REMOVE NA variables and SVL
#pairs(data)
Region <- morpho[,1] #isolate groupings
z <- lda(morpho$clade ~ ., data=data)
pz <- predict(z)
table <- table(pz$class, morpho$clade)
write.csv(table, "dfa_table.csv")
accuracy <- pz$class == morpho$clade
write.csv(accuracy, "accuracy.csv")
pz$class[accuracy==F]
morpho$clade[accuracy==F]
data.frame(original=morpho$clade[accuracy==F], prediction=pz$class[accuracy==F])

scores <-data.frame(morpho$clade, pz$x[])
scores

ldahist(data = pz$x[,1], g=morpho$Region)
plot(pz$x, col = cols, pch=20, cex = 3)
legend(x=2.25,y=4,legend=c(levels(morpho$Region)), col=levels(popcols),pch=c(20,20,20))

pz$x

##to make plot using ggplot
cbPalette3 <- c("#999999", "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")
myCol <- c("#009E73", "#56B4E9", "#E69F00")

dfa <- ggplot(scores, aes(x=LD1, y=LD2, color = scores$morpho.clade)) + theme_bw() + theme(legend.key = element_blank()) +
  scale_color_manual(values = myCol) + stat_ellipse() +
  scale_shape_manual(values=c(20,20,20)) + 
  geom_point(aes(color = scores$morpho.clade), size=5) + 
  theme(legend.title=element_blank(), axis.text=element_text(size=16), 
  axis.title=element_text(size=16), legend.text = element_text(size=20)) + 
  guides(colour = guide_legend(override.aes = list(size=6)))

dfa

ggsave(plot=dfa,height=10,width=10,dpi=300, filename="3regionDFA.pdf", useDingbats=FALSE)
