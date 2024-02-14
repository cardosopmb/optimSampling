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

par(mfrow=c(2,4), cex = 1, lwd = 3, mar = c(2,4,0.2,1))

###################Alpha TD#########################
res = ResAlphaTD[-1,]

plot(res[,11], ylab = 'alpha TD', xlab = 'Accumulation', cex = 0, xlim=c(0,320), ylim = c(0,1))
lines(RndAlphaTDquant[1,], col = "antiquewhite3", lwd = 1)
lines(RndAlphaTDquant[2,], col = "antiquewhite3", lwd = 1)
lines(RndAlphaTDquant[3,], col = "antiquewhite3", lwd = 1)
lines(res[,11], lwd = 1)

#colcodes = c("cadetblue3", "cadetblue4", "chartreuse3", "chartreuse4", "chocolate3", "chocolate4", "antiquewhite4", "brown3", "brown4")
#plot(ResAlphaTD[,7], main = "Samples/method", ylim = c(0,64), ylab = '', cex = 0)
#for (i in c(1:9)){
#	lines(ResAlphaTD[,i], col = colcodes[i], xlab = NULL)
#}
#plot(res[,7], main = "% Samples/method", ylim = c(0,1), ylab = '', cex = 0)
#for (i in c(1:9)){
#	lines(res[,i], col = colcodes[i], ylab = NA)
#}

###################Alpha PD##########################
res = ResAlphaPD[-1,]

plot(res[,11], ylab = 'alpha PD', xlab = '', cex = 0, xlim=c(0,320), ylim = c(0,1))
lines(RndAlphaPDquant[1,], col = "antiquewhite3", lwd = 1)
lines(RndAlphaPDquant[2,], col = "antiquewhite3", lwd = 1)
lines(RndAlphaPDquant[3,], col = "antiquewhite3", lwd = 1)
lines(res[,11], lwd = 1)

#plot(ResAlphaPD[,7], ylim = c(0,64), ylab = '', xlab = 'Samples', cex = 0)
#for (i in c(1:9)){
#  lines(ResAlphaPD[,i], col = colcodes[i])
#}
#plot(res[,7], ylim = c(0,1), ylab = '', xlab = 'Samples', cex = 0)
#for (i in c(1:9)){
#  lines(res[,i], col = colcodes[i])
#}

###################Alpha FD##########################
res = ResAlphaFD[-1,]

plot(res[,11], ylab = 'alpha FD', xlab = 'Samples', cex = 0, xlim=c(0,320), ylim = c(0,1))
lines(RndAlphaFDquant[1,], col = "antiquewhite3", lwd = 1)
lines(RndAlphaFDquant[2,], col = "antiquewhite3", lwd = 1)
lines(RndAlphaFDquant[3,], col = "antiquewhite3", lwd = 1)
lines(res[,11], lwd = 1)

#plot(ResAlphaFD[,7], ylim = c(0,64), ylab = '', xlab = 'Samples', cex = 0)
#for (i in c(1:9)){
#  lines(ResAlphaFD[,i], col = colcodes[i])
#}
#plot(res[,7], ylim = c(0,1), ylab = '', xlab = 'Samples', cex = 0)
#for (i in c(1:9)){
#  lines(res[,i], col = colcodes[i])
#}

###################Alpha Compl (TD Geres & Arrabida, start = 12 PF)##########################
res = ResAlphaCompl

plot(res[,11], ylab = 'alpha TD (constrained)', xlab = 'Samples', cex = 0, ylim = c(0,1))
lines(RndAlphaTDquant[1,-c(1:11)], col = "antiquewhite3", lwd = 1)
lines(RndAlphaTDquant[2,-c(1:11)], col = "antiquewhite3", lwd = 1)
lines(RndAlphaTDquant[3,-c(1:11)], col = "antiquewhite3", lwd = 1)
lines(res[,11], lwd = 1)

#plot(ResAlphaCompl[,7], ylim = c(0,64), ylab = '', xlab = 'Samples', cex = 0)
#for (i in c(1:9)){
#  lines(ResAlphaCompl[,i], col = colcodes[i])
#}
#plot(res[,7], ylim = c(0,1), ylab = '', xlab = 'Samples', cex = 0)
#for (i in c(1:9)){
#  lines(res[,i], col = colcodes[i])
#}

###################Beta TD#########################

res = resBetaTD[-1,]

plot(res[,8], ylab = 'beta TD', xlab = '', cex = 0, xlim=c(0,24), ylim = c(0,1))
lines(RndBetaTDquant[1,], col = "antiquewhite3", lwd = 1)
lines(RndBetaTDquant[2,], col = "antiquewhite3", lwd = 1)
lines(RndBetaTDquant[3,], col = "antiquewhite3", lwd = 1)
lines(res[,8], lwd = 1)

###################Beta PD##########################
res = resBetaPD[-1,]

plot(res[,8], ylab = 'beta PD', xlab = 'Samples', cex = 0, xlim=c(0,24), ylim = c(0,1))
lines(RndBetaPDquant[1,], col = "antiquewhite3", lwd = 1)
lines(RndBetaPDquant[2,], col = "antiquewhite3", lwd = 1)
lines(RndBetaPDquant[3,], col = "antiquewhite3", lwd = 1)
lines(res[,8], lwd = 1)

###################Beta FD##########################
res = resBetaFD[-1,]

plot(res[-1,8], ylab = 'beta FD', xlab = 'Samples', cex = 0, xlim=c(0,24), ylim = c(0,1))
lines(RndBetaFDquant[1,], col = "antiquewhite3", lwd = 1)
lines(RndBetaFDquant[2,], col = "antiquewhite3", lwd = 1)
lines(RndBetaFDquant[3,], col = "antiquewhite3", lwd = 1)
lines(res[,8], lwd = 1)

###################Beta Compl (TD Geres & Arrabida, start = 4 PF)##########################
res = resBetaCompl

plot(c(4:24), res[,8], ylab = 'beta TD (constrained)', xlab = 'Samples', cex = 0, xlim=c(0,24), ylim = c(0,1))
lines(c(4:24), RndBetaTDquant[1,-c(1:3)], col = "antiquewhite3", lwd = 1)
lines(c(4:24), RndBetaTDquant[2,-c(1:3)], col = "antiquewhite3", lwd = 1)
lines(c(4:24), RndBetaTDquant[3,-c(1:3)], col = "antiquewhite3", lwd = 1)
lines(c(4:24), res[,8], lwd = 1)

