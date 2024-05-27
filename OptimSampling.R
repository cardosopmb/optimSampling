#####################################################
##############SPIDER ANALYSES########################
#####################################################

require(abind)
require(ape)
require(BAT)
data(geres)
data(arrabida)
data(guadiana)
data(phylotree)
data(functree)

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
Color = factor(rep(c(rep("grey", 963), rep("darkblue", 24), "darkorange", rep("darkblue", 296)), 4))
TDiv = c(RndAlphaTDquant[1,], RndAlphaTDquant[2,], RndAlphaTDquant[3,], ResAlphaTD[,11])
PDiv = c(RndAlphaPDquant[1,], RndAlphaPDquant[2,], RndAlphaPDquant[3,], ResAlphaPD[,11])
FDiv = c(RndAlphaFDquant[1,], RndAlphaFDquant[2,], RndAlphaFDquant[3,], ResAlphaFD[,11])
CDiv = c(RndAlphaTDquant[1,], RndAlphaTDquant[2,], RndAlphaTDquant[3,], ResAlphaCompl[,11])
Diversity = c(TDiv, PDiv, FDiv, CDiv)

alphaSpiders = data.frame(Dimension, Cost, Color, Diversity)

plot.alphaSpiders = ggplot(data = alphaSpiders) +
  geom_point(data = subset(alphaSpiders, Color == "grey"), aes(x = Cost, y = Diversity, color = Color)) +
  geom_point(data = subset(alphaSpiders, Color == "darkblue"), aes(x = Cost, y = Diversity, color = Color)) +
  geom_point(data = subset(alphaSpiders, Color == "darkorange"), aes(x = Cost, y = Diversity, color = Color, size = 1)) +
  scale_colour_manual(values = c("grey", "darkblue", "darkorange")) +
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
Color = factor(rep(c(rep("grey", 75), rep("darkblue", 6), "darkorange", rep("darkblue", 18)), 4))
TDiv = c(RndBetaTDquant[1,], RndBetaTDquant[2,], RndBetaTDquant[3,], resBetaTD[,8])
PDiv = c(RndBetaPDquant[1,], RndBetaPDquant[2,], RndBetaPDquant[3,], resBetaPD[,8])
FDiv = c(RndBetaFDquant[1,], RndBetaFDquant[2,], RndBetaFDquant[3,], resBetaFD[,8])
CDiv = c(RndBetaTDquant[1,], RndBetaTDquant[2,], RndBetaTDquant[3,], resBetaCompl[,8])
Diversity = c(TDiv, PDiv, FDiv, CDiv)

betaSpiders = data.frame(Dimension, Cost, Color, Diversity)

plot.betaSpiders = ggplot(data = betaSpiders) +
  geom_point(data = subset(betaSpiders, Color == "grey"), aes(x = Cost, y = Diversity, color = Color)) +
  geom_point(data = subset(betaSpiders, Color == "darkblue"), aes(x = Cost, y = Diversity, color = Color)) +
  geom_point(data = subset(betaSpiders, Color == "darkorange"), aes(x = Cost, y = Diversity, color = Color, size = 1)) +
  scale_colour_manual(values = c("grey", "darkblue", "darkorange")) +
  theme(legend.position = "none") +
  facet_wrap(~factor(Dimension, levels = c("Taxonomic", "Phylogenetic", "Functional", "Constrained")), scales = "free", nrow = 1) +
  theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11)) +
  theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) +
  theme(strip.text.x = element_text(size = 14), strip.text.y = element_text(size = 14)) +
  ylab("1 - Bias")
plot.betaSpiders

ggsave("Fig1b.tiff", plot = plot.betaSpiders)
ggsave("Fig1b.jpeg", plot = plot.betaSpiders)



#####################################################
##############MAMMAL ANALYSES########################
#####################################################

#loading the dataset
bats.optim <- read.csv("Bats_optim.csv",h=T)
bats.optim1 <- bats.optim[,3:ncol(bats.optim)]#removing the two first columns

#alpha-sampling
#creating the dataframe with the necessary information regarding sample size and costs
methods.alpha <- data.frame(method = c("mistnets", "audiomoths"),
                       nSamples = c(12, 12),
                       fixCost = c(520, 883.2),
                       varCost = c(384, 1062))

#running the analysis for alpha-sampling
res.bats <- optim.alpha(bats.optim1, methods = methods.alpha)
plot(res.bats$cost, res.bats$diversity)#creating a plot to visualize the benefit-cost ratio
identify(res.bats$cost, res.bats$diversity)#to identify the point with better benefit-cost ratio

#beta-sampling
#isolating the sites to create an array
s1 <- as.matrix(bats.optim[bats.optim$site=="Site_1", 3:ncol(bats.optim)])
s2 <- as.matrix(bats.optim[bats.optim$site=="Site_2", 3:ncol(bats.optim)])
s3 <- as.matrix(bats.optim[bats.optim$site=="Site_3", 3:ncol(bats.optim)])
s4 <- as.matrix(bats.optim[bats.optim$site=="Site_4", 3:ncol(bats.optim)])
s5 <- as.matrix(bats.optim[bats.optim$site=="Site_5", 3:ncol(bats.optim)])
s6 <- as.matrix(bats.optim[bats.optim$site=="Site_6", 3:ncol(bats.optim)])
s7 <- as.matrix(bats.optim[bats.optim$site=="Site_7", 3:ncol(bats.optim)])
s8 <- as.matrix(bats.optim[bats.optim$site=="Site_8", 3:ncol(bats.optim)])
s9 <- as.matrix(bats.optim[bats.optim$site=="Site_9", 3:ncol(bats.optim)])
s10 <- as.matrix(bats.optim[bats.optim$site=="Site_10", 3:ncol(bats.optim)])
s11 <- as.matrix(bats.optim[bats.optim$site=="Site_11", 3:ncol(bats.optim)])
s12 <- as.matrix(bats.optim[bats.optim$site=="Site_12", 3:ncol(bats.optim)])

#creating the array
bats.beta <- array(c(s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12), c(2, ncol(bats.optim)-2, 12))
colnames(bats.beta) <- colnames(bats.optim)[3:ncol(bats.optim)]

#creating the dataframe with the necessary information regarding sample size and costs
methods.beta <- data.frame(method = c("mistnets","audiomoths"),
                            nSamples = c(1, 1),
                            fixCost = c(520, 883.2),
                            varCost = c(384, 1062))

#running the analysis for beta-sampling
res.bats.beta <- optim.beta(bats.beta, methods = methods.beta)
plot(res.bats.beta$cost, res.bats.beta$diversity)#creating a plot to visualize the benefit-cost ratio
identify(res.bats.beta$cost, res.bats.beta$diversity)#to identify the point with better benefit-cost ratio

#phylogenetic analysis
#loading the dataset
bats.phylo <- read.csv("bats_phylo.csv", sep=";", header=T, row.names = 1)
tree <- hclust(dist(bats.phylo, method="binary"), method="average")

#running the analysis for alpha-sampling
res.phylo.bat <- optim.alpha(bats.optim1,tree, methods = methods.alpha)
plot(res.phylo.bat$cost, res.phylo.bat$diversity)#creating a plot to visualize the benefit-cost ratio
identify(res.phylo.bat$cost, res.phylo.bat$diversity)#to identify the point with better benefit-cost ratio

#running the analysis for beta-sampling
phylo.bats.beta=optim.beta(bats.beta, tree, methods = methods.beta)
plot(phylo.bats.beta$cost, phylo.bats.beta$diversity)#creating a plot to visualize the benefit-cost ratio
identify(phylo.bats.beta$cost, phylo.bats.beta$diversity)#to identify the point with better benefit-cost ratio

#functional analysis
#loading the dataset
bats.func <- read.csv("Bats_optim_func.csv", h=T, sep=";", dec=",", row.names = 1)
func.tree.bats <- tree.build(gawdis(bats.func))

#running the analysis for alpha-sampling
res.bats.func <- optim.alpha(bats.optim1, func.tree.bats, methods = methods.alpha)
plot(res.bats.func$cost, res.bats.func$diversity)#creating a plot to visualize the benefit-cost ratio
identify(res.bats.func$cost,res.bats.func$diversity)#to identify the points with better benefit-cost ratio

#running the analysis for beta-sampling
func.bats.beta <- optim.beta(bats.beta, func.tree.bats, methods = methods1.beta)
plot(func.bats.beta$cost, func.bats.beta$diversity)#creating a plot to visualize the benefit-cost ratio
identify(func.bats.beta$cost, func.bats.beta$diversity)#to identify the point with better benefit-cost ratio


#making the plots
library(ggplot2)
# rescaling the diversity values to obtain values from 0 to 1
res.bats1 <- res.bats; res.bats1$diversity <- res.bats$diversity/max(res.bats$diversity)#TD
res.bats.func1 <- res.bats.func; res.bats.func1$diversity <- res.bats.func$diversity/max(res.bats.func$diversity)#PD
res.phylo.bat1 <- res.phylo.bat; res.phylo.bat1$diversity <- res.phylo.bat$diversity/max(res.phylo.bat$diversity)#FD

#creating a dataframe with the rescaled results
alpha=rbind(res.bats1, res.phylo.bat1, res.bats.func1)

#creating a vector called 'best', which identifies the optimal choices assessed through the plots above
n <- nrow(res.bats1)
best <- rep("no", nrow(alpha)); best[c(128,
                                         n+46,n+56, n+60, n+65, n+71, n+77, n+80, n+86 ,n+89, n+91, n+95, 
                                         2*n+54,2*n+65,2*n+75,2*n+83,2*n+92,2*n+98,2*n+105,2*n+112,2*n+117)]="yes"
#note that there are several optimal choices for bat PD
alpha$best <- best
colnames(alpha) <- c("Method1","Method2","Cost","Diversity","Color")

#creating columns to identify the results
alpha$type=rep(c('1Taxonomic', '2Phylogenetic', '3Functional'),
               c(nrow(res.bats), nrow(res.phylo.bat), nrow(res.bats.func)))

#plot
#the same as above, but for beta diversity
beta <- rbind(res.bats.beta, phylo.bats.beta, func.bats.beta)
best.beta <- rep("no",nrow(beta)); best.beta[c(3, 7, 11)] <- "yes"

beta$best <- best.beta
colnames(beta) <- c("Method1", "Method2", "Cost", "Diversity", "Color")

beta$type=rep(c('1Taxonomic', '2Phylogenetic', '3Functional'),
               c(nrow(res.bats.beta), nrow(phylo.bats.beta), nrow(func.bats.beta)))

#combining the dataframes and including a column to specify the type of sampling
all <- rbind(alpha,beta)
all$sampling <- rep(c("α-sampling","β-sampling"),c(nrow(alpha),nrow(beta)))

#specifying the lables
type.labels=c(`1Taxonomic`='Taxonomic',
              `2Phylogenetic`="Phylogenetic",
              `3Functional`='Functional',
              `α-sampling`="Inventorying",
              `β-sampling`="Monitoring")

#creating the plot
plot=ggplot(data = all) + geom_point(aes(x = Cost, y = Diversity, color = Color, size = Color))+
  scale_color_manual(values = c("darkblue", "darkorange"), labels = NULL) +
  scale_size_manual(values=c(1, 2)) +
  theme(legend.position = "none") + ylab("Average proportion of diversity (α) or 1 - bias (β)") +
  facet_grid(vars(type), vars(sampling), scales = "free", labeller = as_labeller(type.labels)) +
  theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11)) +
  theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14))+
  theme(strip.text.x = element_text(size = 14), strip.text.y = element_text(size = 14))

#saving the plot
ggsave("plot.optim.tiff",plot = plot)
ggsave("plot.optim.jpg",plot = plot)



