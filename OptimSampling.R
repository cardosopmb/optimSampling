require(abind)
require(ape)
require(BAT)
data(geres)
data(arrabida)
data(guadiana)
data(phylotree)
data(functree)
geres
arrabida
guadiana
phylotree
functree

#####################################################
#####################ANALYSES########################
#####################################################

r = 1000

##methods
Methods320 <- data.frame(Methods = c("AD", "AN", "BD", "BN", "GD", "GN", "PF", "SD", "SN"), Samples = c(32, 32, 32, 32, 32, 32, 64, 32, 32))
Methods24 <- data.frame(Methods = c("AN", "BD", "BN", "PF", "SD", "SN"), Samples = c(4, 2, 2, 12, 2, 2))

##alpha optimal, sequential, cost always the same for every method
dataAlpha <- abind(geres, arrabida, along = 3)
ResAlphaTD <- optim.alpha(dataAlpha, methods = Methods320, seq = TRUE, runs = r)
ResAlphaPD <- optim.alpha(dataAlpha, phylotree, methods = Methods320, seq = TRUE, runs = r)
ResAlphaFD <- optim.alpha(dataAlpha, functree, methods = Methods320, seq = TRUE, runs = r)
ResAlphaCompl <- optim.alpha(dataAlpha, methods = Methods320, base = c(0,0,0,0,0,0,12,0,0), seq = TRUE, runs = r)

##alpha random
RndAlphaTD = matrix(rep(0,320*r),nrow=r)
RndAlphaPD = RndAlphaTD
RndAlphaFD = RndAlphaTD
methods = data.frame(Methods = "A", Samples = 320)
for(i in 1:r){
  cat(i, ",")
  for(j in 1:320){
    RndAlphaTD[i,j] <- as.vector(optim.alpha.stats(dataAlpha, methods = methods, samples = c(j), runs = 1))
    RndAlphaPD[i,j] <- as.vector(optim.alpha.stats(dataAlpha, tree = phylotree,   methods = methods, samples = c(j), runs = 1))
    RndAlphaFD[i,j] <- as.vector(optim.alpha.stats(dataAlpha, tree = functree, methods = methods, samples = c(j), runs = 1))
  }
}
RndAlphaTDquant <- apply(RndAlphaTD, 2, quantile, probs = c(0.025,0.5,0.975))
RndAlphaPDquant <- apply(RndAlphaPD, 2, quantile, probs = c(0.025,0.5,0.975))
RndAlphaFDquant <- apply(RndAlphaFD, 2, quantile, probs = c(0.025,0.5,0.975))

#################################################
#BETA
#################################################

setwd("C:/Users/ungol/OneDrive - Universidade de Lisboa/Ambiente de Trabalho/Papers/OptimSampling")
phylotree = read.tree(file = "PDTree.nex")
functree = read.tree(file = "FDTree.nex") 
Methods24 <- data.frame(Methods = c("AN", "BD", "BN", "PF", "SD", "SN"), Samples = c(4, 2, 2, 12, 2, 2))

##site data (16 from Spain)
dataSites <- read.csv("AbunData.csv", header=TRUE)
colnames(dataSites)[1] = "Aelurillus_luctuosus"

dataBeta <- dataSites[1:24,]
for(i in 2:16){
  dataBeta <- abind(dataBeta, dataSites[((i-1)*24+1):(i*24),], along = 3)
}

##beta optimal
r = 100
resBetaTD <- optim.beta(dataBeta, methods = Methods24, seq = TRUE, runs = r)
resBetaPD <- optim.beta(dataBeta, phylotree, methods = Methods24, seq = TRUE, runs = r)
resBetaFD <- optim.beta(dataBeta, functree, methods = Methods24, seq = TRUE, runs = r)
resBetaCompl <- optim.beta(dataBeta, methods = Methods24, base = c(0,0,0,4,0,0), seq = TRUE, runs = r)

##beta random
RndBetaTD = matrix(rep(0,24*r),nrow=r)
RndBetaPD = RndBetaTD
RndBetaFD = RndBetaTD
methods = data.frame(Methods = "A", Samples = 24)
for(i in 1:r){
  for(j in 1:24){
    RndBetaTD[i,j] <- as.vector(optim.beta.stats(dataBeta, methods = methods, samples = c(j), runs = 1))
    RndBetaPD[i,j] <- as.vector(optim.beta.stats(dataBeta, tree = phylotree, methods = methods, samples = c(j), runs = 1))
    RndBetaFD[i,j] <- as.vector(optim.beta.stats(dataBeta, tree = functree, methods = methods, samples = c(j), runs = 1))
  }
  cat(i, ",")
} 
RndBetaTDquant <- apply(RndBetaTD, 2, quantile, probs = c(0.025,0.5,0.975))
RndBetaPDquant <- apply(RndBetaPD, 2, quantile, probs = c(0.025,0.5,0.975))
RndBetaFDquant <- apply(RndBetaFD, 2, quantile, probs = c(0.025,0.5,0.975))

#####################################################
######################GRAPHS#########################
#####################################################

library(ggplot2)

ResAlphaCompl = rbind(matrix(0, nrow = 12, ncol = ncol(ResAlphaCompl)), ResAlphaCompl)
RndAlphaTDquant = cbind(rep(0, 3), RndAlphaTDquant)
RndAlphaPDquant = cbind(rep(0, 3), RndAlphaPDquant)
RndAlphaFDquant = cbind(rep(0, 3), RndAlphaFDquant)

Dimension = c(rep("Taxonomic", 1284), rep("Phylogenetic", 1284), rep("Functional", 1284), rep("Constrained", 1284))
Cost = rep(0:320, 16)
Color = factor(rep(c(rep("darkblue", 963), rep("darkorange", 321)), 4))
TDiv = c(RndAlphaTDquant[1,], RndAlphaTDquant[2,], RndAlphaTDquant[3,], ResAlphaTD[,11])
PDiv = c(RndAlphaPDquant[1,], RndAlphaPDquant[2,], RndAlphaPDquant[3,], ResAlphaPD[,11])
FDiv = c(RndAlphaFDquant[1,], RndAlphaFDquant[2,], RndAlphaFDquant[3,], ResAlphaFD[,11])
CDiv = c(RndAlphaTDquant[1,], RndAlphaTDquant[2,], RndAlphaTDquant[3,], ResAlphaCompl[,11])
Diversity = c(TDiv, PDiv, FDiv, CDiv)

alphaSpiders = data.frame(Dimension, Cost, Color, Diversity)

plot.alphaSpiders = ggplot(data = alphaSpiders) +
  geom_point(aes(x = Cost, y = Diversity, color = Color)) +
  scale_colour_manual(values = c("darkblue", "darkorange")) +
  theme(legend.position = "none") +
  facet_wrap(~factor(Dimension, levels = c("Taxonomic", "Phylogenetic", "Functional", "Constrained")), scales = "free", nrow = 1) +
  theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11)) +
  theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) +
  theme(strip.text.x = element_text(size = 14), strip.text.y = element_text(size = 14))
plot.alphaSpiders

ggsave("Fig1a.tiff", plot = plot.alphaSpiders)
ggsave("Fig1a.jpeg", plot = plot.alphaSpiders)

###########Beta

resBetaCompl = rbind(matrix(0, nrow = 4, ncol = ncol(resBetaCompl)), resBetaCompl)
RndBetaTDquant = cbind(rep(0, 3), RndBetaTDquant)
RndBetaPDquant = cbind(rep(0, 3), RndBetaPDquant)
RndBetaFDquant = cbind(rep(0, 3), RndBetaFDquant)

Dimension = c(rep("Taxonomic", 100), rep("Phylogenetic", 100), rep("Functional", 100), rep("Constrained", 100))
Cost = rep(0:24, 16)
Color = factor(rep(c(rep("darkblue", 75), rep("darkorange", 25)), 4))
TDiv = c(RndBetaTDquant[1,], RndBetaTDquant[2,], RndBetaTDquant[3,], resBetaTD[,8])
PDiv = c(RndBetaPDquant[1,], RndBetaPDquant[2,], RndBetaPDquant[3,], resBetaPD[,8])
FDiv = c(RndBetaFDquant[1,], RndBetaFDquant[2,], RndBetaFDquant[3,], resBetaFD[,8])
CDiv = c(RndBetaTDquant[1,], RndBetaTDquant[2,], RndBetaTDquant[3,], resBetaCompl[,8])
Diversity = c(TDiv, PDiv, FDiv, CDiv)

betaSpiders = data.frame(Dimension, Cost, Color, Diversity)

plot.betaSpiders = ggplot(data = betaSpiders) +
  geom_point(aes(x = Cost, y = Diversity, color = Color)) +
  scale_colour_manual(values = c("darkblue", "darkorange")) +
  theme(legend.position = "none") +
  facet_wrap(~factor(Dimension, levels = c("Taxonomic", "Phylogenetic", "Functional", "Constrained")), scales = "free", nrow = 1) +
  theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11)) +
  theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) +
  theme(strip.text.x = element_text(size = 14), strip.text.y = element_text(size = 14)) +
  ylab("1 - Bias")
plot.betaSpiders

ggsave("Fig1b.tiff", plot = plot.betaSpiders)
ggsave("Fig1b.jpeg", plot = plot.betaSpiders)
