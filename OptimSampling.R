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
geres
arrabida
guadiana
phylotree
functree

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

library(BAT)
library(gawdis)
library(ape)
library(picante)
dat.optim=read.csv("dados_optim.csv",h=T)
dat.optim1=dat.optim[order(dat.optim$Method),3:30]
methods=data.frame(method=c("camera.trap","search"),
                   nSamples=c(10,10),
                   fixCost=c(192,384),
                   varCost=c(359.51,0))#custo fixo e custo variavel foram inseridos aqui

res=optim.alpha(dat.optim1, methods = methods)
plot(res$cost,res$diversity)#grafico simples para visualizar os resultados
identify(res$cost,res$diversity)#melhor custo-beneficio foi o 99 (so transectos)


#beta
#separando os dados por transecto para criar o array
t1=as.matrix(dat.optim[dat.optim$Transecto=="T1",3:30])
t2=as.matrix(dat.optim[dat.optim$Transecto=="T2",3:30])
t3=as.matrix(dat.optim[dat.optim$Transecto=="T3",3:30])
t4=as.matrix(dat.optim[dat.optim$Transecto=="T4",3:30])
t5=as.matrix(dat.optim[dat.optim$Transecto=="T5",3:30])
t6=as.matrix(dat.optim[dat.optim$Transecto=="T6",3:30])
t7=as.matrix(dat.optim[dat.optim$Transecto=="T7",3:30])
t8=as.matrix(dat.optim[dat.optim$Transecto=="T8",3:30])
t9=as.matrix(dat.optim[dat.optim$Transecto=="T9",3:30])
t10=as.matrix(dat.optim[dat.optim$Transecto=="T10",3:30])

dat.beta=array(c(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10),c(2,28,10))
colnames(dat.beta)=colnames(dat.optim)[3:30]

methods.beta=data.frame(method=c("camera.trap","search"),
                        nSamples=c(1,1),
                        fixCost=c(192,384),
                        varCost=c(359.51,0))

res.beta=optim.beta(dat.beta,methods=methods.beta)
plot(res.beta$cost,res.beta$diversity)#grafico simples para visualizar os resultados

#importacao de tracos funcionaist1=dat.optim[dat.optim$Transecto=="T1",]
func.data=read.csv2("functional traits.csv",h=T,row.names = 1)
#combinando variaveis de dieta (presas e partes reprodutivas de plantas)
func.data$Diet.prey=rowSums(func.data[,1:6])
func.data$Diet.rep=rowSums(func.data[,7:9])
func.data=func.data[,-c(1:9)]#removendo colunas excedentes
func.data$BodyMass.Value=log(as.numeric(func.data$BodyMass.Value))

#criando arevore funcional com gawdis e tree.build
func.tree <- tree.build(gawdis(func.data,groups=c(1,2,3,3,3,4,1,1)))
res.func=optim.alpha(dat.optim1,func.tree, methods = methods)
plot(res.func$cost,res.func$diversity)
identify(res.func$cost,res.func$diversity)

#beta
res.beta.func=optim.beta(dat.beta,func.tree,methods=methods.beta)
plot(res.beta.func$cost,res.beta.func$diversity)

#analise filogenetica
phylotree=consensus(read.nexus("output.nex"),rooted = T)#importando a arvore filogenetica e obtendo uma  uma arvore de consenso a partir das 1000 arvores
plot(phylotree)
length(phylotree$tip.label);ncol(dat.optim2)
dat.optim2=dat.optim1#criando um novo objeto para ajustar os nomes das especies
#isso foi necessario pq o nome das especies nos dados e na arvore filogenetica tem algumas diferencas
colnames(dat.optim2)=c(colnames(dat.optim1)[1:5],"Leopardus_wiedii",colnames(dat.optim1)[7:12],
                       "Sciurus_aestuans",colnames(dat.optim1)[14:22],
                       "Hydrochoerus_hydrochaeris",colnames(dat.optim1)[24:27],"Canis_lupus")
#nao consegui inserir a arvore de consenso na analise, mas usei a arvore 6452, que tem a mesma topografia
res.phylo=optim.alpha(dat.optim2,phylotreex$tree_6452, methods = methods)
plot(res.phylo$cost,res.phylo$diversity)
identify(res.phylo$cost,res.phylo$diversity)


#beta
#modificando o nome das especies
colnames(dat.beta)=c(colnames(dat.optim1)[1:5],"Leopardus_wiedii",colnames(dat.optim1)[7:12],
                     "Sciurus_aestuans",colnames(dat.optim1)[14:22],
                     "Hydrochoerus_hydrochaeris",colnames(dat.optim1)[24:27],"Canis_lupus")
res.beta.phylo=optim.beta(dat.beta,phylotreex$tree_6452,methods=methods.beta)
plot(res.beta.phylo$cost,res.beta.phylo$diversity)


########################
# análise com morcegos #
########################
bats.optim=read.csv("Morcegos_optim.csv",h=T)
bats.optim1=bats.optim[,3:ncol(bats.optim)]

methods1=data.frame(method=c("mistnets","audiomoths"),
                    nSamples=c(12,12),
                    fixCost=c(520,883.2),
                    varCost=c(384,1062))

res.bats=optim.alpha(bats.optim1, methods = methods1)
plot(res.bats$cost,res.bats$diversity)
identify(res.bats$cost,res.bats$diversity)

#separando os sites para fazer o array (div beta)
s1=as.matrix(bats.optim[bats.optim$site=="Site_1",3:ncol(bats.optim)])
s2=as.matrix(bats.optim[bats.optim$site=="Site_2",3:ncol(bats.optim)])
s3=as.matrix(bats.optim[bats.optim$site=="Site_3",3:ncol(bats.optim)])
s4=as.matrix(bats.optim[bats.optim$site=="Site_4",3:ncol(bats.optim)])
s5=as.matrix(bats.optim[bats.optim$site=="Site_5",3:ncol(bats.optim)])
s6=as.matrix(bats.optim[bats.optim$site=="Site_6",3:ncol(bats.optim)])
s7=as.matrix(bats.optim[bats.optim$site=="Site_7",3:ncol(bats.optim)])
s8=as.matrix(bats.optim[bats.optim$site=="Site_8",3:ncol(bats.optim)])
s9=as.matrix(bats.optim[bats.optim$site=="Site_9",3:ncol(bats.optim)])
s10=as.matrix(bats.optim[bats.optim$site=="Site_10",3:ncol(bats.optim)])
s11=as.matrix(bats.optim[bats.optim$site=="Site_11",3:ncol(bats.optim)])
s12=as.matrix(bats.optim[bats.optim$site=="Site_12",3:ncol(bats.optim)])

bats.beta=array(c(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12),c(2,ncol(bats.optim)-2,12))
colnames(bats.beta)=colnames(bats.optim)[3:ncol(bats.optim)]

methods1.beta=data.frame(method=c("mistnets","audiomoths"),
                         nSamples=c(1,1),
                         fixCost=c(520,883.2),
                         varCost=c(384,1062))

res.bats.beta=optim.beta(bats.beta, methods = methods1.beta)
plot(res.bats.beta$cost,res.bats.beta$diversity)

#funcional
bats.func=read.csv("Morcegos_optim_func.csv",h=T,dec=",",row.names = 1)
func.tree.bats <- tree.build(gawdis(bats.func))

res.bats.func=optim.alpha(bats.optim1,func.tree.bats, methods = methods1)
plot(res.bats.func$cost,res.bats.func$diversity)
identify(res.bats.func$cost,res.bats.func$diversity)#91, 101, 108, 118

func.bats.beta=optim.beta(bats.beta, func.tree.bats,methods = methods1.beta)
plot(func.bats.beta$cost,func.bats.beta$diversity)

#filogenetica
#aqui usei uma arvore fornecida pelo William com o pacote 'picante'
bats_super_phylo <- read.csv("bats_phylo_Renato.csv", sep=";", header=T, row.names = 1)
tree <- hclust(dist(bats_super_phylo, method="binary"), method="average")
plot(tree)

res.phylo.bat=optim.alpha(bats.optim1,tree, methods = methods1)
plot(res.phylo.bat$cost,res.phylo.bat$diversity)
identify(res.phylo.bat$cost,res.phylo.bat$diversity)

phylo.bats.beta=optim.beta(bats.beta, tree,methods = methods1.beta)
plot(phylo.bats.beta$cost,phylo.bats.beta$diversity)


##################
#### graficos ####
##################

library(ggplot2)

#criando novos objetos com os resultados para poder padronizar os valores de diversidade entre 0 e 1
res1=res;res1$diversity=res$diversity/max(res$diversity)#TD de mamiferos
res.func1=res.func;res.func1$diversity=res.func$diversity/max(res.func$diversity)#FD de mamiferos
res.phylo1=res.phylo;res.phylo1$diversity=res.phylo$diversity/max(res.phylo$diversity)#PD de mamiferos
res.bats1=res.bats;res.bats1$diversity=res.bats$diversity/max(res.bats$diversity)#TD de morcegos
res.phylo.bat1=res.phylo.bat;res.phylo.bat1$diversity=res.phylo.bat$diversity/max(res.phylo.bat$diversity)#FD de morcegos
res.bats.func1=res.bats.func;res.bats.func1$diversity=res.bats.func$diversity/max(res.bats.func$diversity)#PD de morcegos

#mamiferos
alfa=rbind(res1,res.func1,res.phylo1)
best=rep("black",nrow(alfa));best[c(99,nrow(res)+99,(2*nrow(res))+101)]="red"
#best e um vetor de cores que codifica quais sao as melhores combinacoes de metodos. Foi visto a partir dos graficos simples criados mais acima, usando a funcao 'identify'
alfa$best=best

#morcegos
alfa1=rbind(res.bats1,res.bats.func1,res.phylo.bat1)
n=nrow(res.bats1)
best1=rep("black",nrow(alfa1));best1[c(127,
                                       n+118,
                                       2*n+56,2*n+60,2*n+64,2*n+72,2*n+77,2*n+80,2*n+83,2*n+84,
                                       2*n+94,2*n+98,2*n+100,2*n+102,2*n+106,2*n+107,
                                       2*n+92,2*n+89,2*n+115,2*n+119)]="red"
#no caso de morcegos, houve varias opcoes com boa relacao custo-beneficio
alfa1$best=best1
colnames(alfa1)=colnames(alfa)
alfa2=rbind(alfa,alfa1)#juntando os resultados de mamiferos e morcegos para criar o conjunto de dados
colnames(alfa2)=c("Method1","Method2","Cost","Diversity","Color")

#criando colunas de identificao dos dados
alfa2$taxon=rep(c('Mammals','Bats'),c(nrow(alfa),nrow(alfa1)))
alfa2$type=rep(c('1Taxonomic','2Functional',"3Phylogenetic",'1Taxonomic','2Functional',"3Phylogenetic"),
               c(nrow(res),nrow(res.func),nrow(res.phylo),nrow(res.bats),nrow(res.bats.func),nrow(res.phylo.bat)))

#isso aqui foi feito para ajustar a ordem em que as informacoes aparecem no ggplot (senao, fica na ordem alfabetica)
type.labels=c(`1Taxonomic`='Taxonomic',
              `2Functional`='Functional',
              `3Phylogenetic`="Phylogenetic",
              `Mammals`="Mammals",
              `Bats`="Bats")


plot.alpha=ggplot(data = alfa2)+geom_point(aes(x=Cost,y=Diversity,color=Color))+
  scale_color_manual(values=c("darkblue", "darkorange"),labels=NULL)+theme(legend.position = "none")+
  facet_grid(vars(taxon), vars(type), scales="free",labeller=as_labeller(type.labels))+
  theme(axis.text.x=element_text(size=11),axis.text.y=element_text(size=11))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  theme(strip.text.x=element_text(size=14),strip.text.y=element_text(size=14))
plot.alpha



beta1=rbind(res.beta,res.beta.func,res.beta.phylo)
best.beta=rep("darkblue",nrow(beta1))
best.beta[c(3,7,11)]="orange"
beta2=rbind(res.bats.beta,func.bats.beta,phylo.bats.beta)
best.beta2=rep("darkblue",nrow(beta2))
best.beta2[c(3,7,11)]="orange"
colnames(beta2)=colnames(beta1)
beta1$best=best.beta
beta2$best=best.beta2
beta3=rbind(beta1,beta2)
colnames(beta3)=c("Method1","Method2","Cost","Diversity","Color")

beta3$taxon=rep(c('Mammals','Bats'),c(nrow(beta1),nrow(beta2)))
beta3$type=rep(c('1Taxonomic','2Functional',"3Phylogenetic",'1Taxonomic','2Functional',"3Phylogenetic"),
               c(nrow(res.beta),nrow(res.beta.func),nrow(res.beta.phylo),
                 nrow(res.bats.beta),nrow(func.bats.beta),nrow(phylo.bats.beta)))

plot.beta=ggplot(data = beta3)+geom_point(aes(x=Cost,y=Diversity,color=Color))+
  scale_color_manual(values=c("darkblue", "darkorange"),labels=NULL)+
  theme(legend.position = "none")+ylab("1-Bias")+
  facet_grid(vars(taxon), vars(type), scales="free",labeller =as_labeller(type.labels) )+
  theme(axis.text.x=element_text(size=11),axis.text.y=element_text(size=11))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  theme(strip.text.x=element_text(size=14),strip.text.y=element_text(size=14))
plot.beta

ggsave("plot.alpha.tiff",plot = plot.alpha)
ggsave("plot.beta.tiff",plot = plot.beta)
ggsave("plot.alpha.jpg",plot = plot.alpha)
ggsave("plot.beta.jpg",plot = plot.beta)


