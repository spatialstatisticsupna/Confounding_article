rm(list=ls())

library(INLA)
library(RColorBrewer)
library(rgdal)
library(spdep)
library(tmap)


#########################################
## Read the dowry death mortality data ##
#########################################
data <- read.table(file="../data/DowryDeaths_UttarPradesh.txt", header=TRUE)
str(data)

S <- length(unique(data$dist))
T <- length(unique(data$year))

t.from <- min(data$year)
t.to <- max(data$year)

## Order the data by year and area ##
data <- data[order(data$year,data$dist),]
head(data)

## Expected number of cases and standardized mortality ratio (SMR)
data$E <- data$pop_linear*(sum(data$obs)/sum(data$pop_linear))
data$SMR <- data$obs/data$E
# all.equal(sum(data$obs),sum(data$E))

Data <- data.frame(O=data$obs, E=data$E, SMR=data$SMR,
                   X1=scale(data$x1), # Sex ratio (nº of woman per 100,000 men)
                   X2=scale(data$x2), # Population density (nº of people per square kilometer)
                   X3=scale(data$x3), # Female literacy rate (%)
                   X4=scale(data$x4), # Per capita income (PCI)
                   X5=scale(data$x5), # Murder (per 100,000 people)
                   X6=scale(data$x6), # Burglary (per 100,000 people)
                   ID.area=rep(1:S,T), ID.year=rep(1:T,each=S), ID.area.year=seq(1,S*T))

Beta.df <- as.matrix(Data[,paste("X",1:6,sep="")])


##################################################
## Define spatial and temporal structure matrix ##
##################################################
g <- inla.read.graph("../data/Uttar_Pradesh_nb.graph")
Qs <- matrix(0, g$n, g$n)
for (i in 1:g$n){
  Qs[i,i]=g$nnbs[[i]]
  Qs[i,g$nbs[[i]]]=-1
}

Dm <- diff(diag(T),differences=1)
Qt <- t(Dm)%*%Dm

Qst <- kronecker(Qt,Qs)


###############################
## Read the cartography file ##
###############################
carto_up <- readOGR("../data/carto_up")
plot(carto_up,axes=T)

## Add the neighborhood graph ##
W <- -Qs
diag(W) <- 0
carto.nb <- mat2listw(W)$neighbours
plot(carto.nb, coordinates(carto_up), pch=19, cex=0.8, col="blue", add=TRUE)


#####################################################################################
## Figure: matplot with the temporal evolution of the covariates for each district ##
#####################################################################################
par(mfrow=c(2,3), pty="s")

var.name <- c("Sex ratio","Population density","Female literacy rate","Per capita income","Murder rate","Burglary rate")

for(i in 1:6){
  matplot(t(matrix(Beta.df[,i],S,T,byrow=F)), type="l", main=var.name[i], ylab="", xaxt="n")
  axis(1, at=seq(1,T,2), labels=seq(t.from,t.to,2))
}
par(mfrow=c(1,1))


##########################################################################################
## Load the models fitted with INLA and PQL, which can be downloaded from:              ##
##  - https://emi-sstcdapp.unavarra.es/Confounding_article/data/DataAnalysis_INLA.Rdata ##
##  - https://emi-sstcdapp.unavarra.es/Confounding_article/data/DataAnalysis_PQL.Rdata  ##
##########################################################################################
load("DataAnalysis_INLA.Rdata")
Models.INLA <- list(Model1,Model2,Model3b,Model4)

load("DataAnalysis_PQL.Rdata")
Models.PQL <- list(Model1,Model2,Model3b,Model4)

if(!file.exists("figures")) {
  dir.create(file.path(getwd(), "figures"))
}


#######################################################################################
## Table 1: Posterior mean, posterior standard deviation, and 95% credible intervals ##
##          of the fixed effects for models fitted with INLA, and point estimates,   ##
##          standard errors and 95% confidence intervals obtained with PQL.          ##
#######################################################################################
p <- 7 ## number of fixed effects (including the intercept)

## INLA models ##
INLA.fixed <- lapply(Models.INLA, function(x) x$summary.fixed[,c(1,2,3,5)])

Table1.INLA <- vector("list",p)
names(Table1.INLA) <- c("Intercept",var.name)

for(i in 1:p){
  Table1.INLA[[i]] <- round(do.call(rbind, lapply(INLA.fixed, function(x) x[i,])),4)
  rownames(Table1.INLA[[i]]) <- c("ST1","ST2","ST3","ST4")
}

print(Table1.INLA[-1])


## PQL models ##
PQL.fixed <- list(data.frame(mean=Models.PQL[[1]]$coefficients,
                             sd=sqrt(diag(vcov(Models.PQL[[1]]))),
                             q.025=confint(Models.PQL[[1]])[,1],
                             q.975=confint(Models.PQL[[1]])[,2]))

PQL.fixed <- append(PQL.fixed, lapply(Models.PQL[-1], function(x) data.frame(x$param[1:p,],
                                                                             q.025=x$param[1:p,1]-qnorm(0.975)*x$param[1:p,2],
                                                                             q.975=x$param[1:p,1]+qnorm(0.975)*x$param[1:p,2])))
Table1.PQL <- vector("list",p)
names(Table1.PQL) <- c("Intercept",var.name)

for(i in 1:p){
  Table1.PQL[[i]] <- round(do.call(rbind, lapply(PQL.fixed, function(x) x[i,])),4)
  rownames(Table1.PQL[[i]]) <- c("ST1","ST2","ST3","ST4")
}

print(Table1.PQL[-1])


####################################################################
## Table 2: Model selection criteria and computing time (seconds) ## 
##          for spatio-temporal models fit with INLA and PQL      ##
####################################################################

## INLA model ##
Table2.INLA <- data.frame(mean.deviance=round(unlist(lapply(Models.INLA, function(x) x$dic$mean.deviance)),2),
                          p.eff=round(unlist(lapply(Models.INLA, function(x) x$dic$p.eff)),2),
                          DIC=round(unlist(lapply(Models.INLA, function(x) x$dic$dic)),2),
                          Time=round(unlist(lapply(Models.INLA, function(x) x$cpu.used[4]))),
                          row.names=c("Model ST1","Model ST2","Model ST3","Model ST4"))
print(Table2.INLA)


## PQL model ##
Table2.PQL <- data.frame(Deviance=round(Models.PQL[[1]]$deviance-2*sum(log(dpois(Data$O,Data$O))),2),
                         Df=length(Data$O)-Models.PQL[[1]]$df.residual,
                         AIC=round(Models.PQL[[1]]$aic,2),
                         Time=NA)

Table2.PQL <- rbind(Table2.PQL,
                    data.frame(Deviance=round(unlist(lapply(Models.PQL[-1], function(x) x$Satured.deviance+x$Devian)),2),
                               Df=round(unlist(lapply(Models.PQL[-1], function(x) x$Df)),2),
                               AIC=round(unlist(lapply(Models.PQL[-1], function(x) x$Satured.deviance+x$AIC)),2),
                               Time=round(unlist(lapply(Models.PQL[-1], function(x) x$cpu.time.elapsed)))))
rownames(Table2.PQL) <- c("Model ST1","Model ST2","Model ST3","Model ST4")

print(Table2.PQL)


###############################################################################################
## Figure 1: Scatter plots of relative risk estimates obtained from Models ST2, ST3 and ST4. ##
###############################################################################################

## INLA models ##
graphics.off()

pdf("figures/Figure1a.pdf", width=12, height=4)
par(mfrow=c(1,3), pty="s")

plot(Models.INLA[[2]]$summary.fitted.values$mean, Models.INLA[[3]]$summary.fitted.values$mean[1:(S*T)],
     main="INLA. Model ST2 vs. Model ST3", cex.main=1.5, xlab="Model ST2", ylab="Model ST3", cex.lab=1.5)
lines(c(0,3),c(0,3))

plot(Models.INLA[[2]]$summary.fitted.values$mean, Models.INLA[[4]]$summary.fitted.values$mean[1:(S*T)],
     main="INLA. Model ST2 vs. Model ST4", cex.main=1.5, xlab="Model ST2", ylab="Model ST4", cex.lab=1.5)
lines(c(0,3),c(0,3))

plot(Models.INLA[[3]]$summary.fitted.values$mean[1:(S*T)], Models.INLA[[4]]$summary.fitted.values$mean[1:(S*T)],
     main="INLA. Model ST3 vs. Model ST4", cex.main=1.5, xlab="Model ST3", ylab="Model ST4", cex.lab=1.5)
lines(c(0,3),c(0,3))

dev.off()


## PQL models ##
graphics.off()

pdf("figures/Figure1b.pdf", width=12, height=4)
par(mfrow=c(1,3), pty="s")

plot(Models.PQL[[2]]$risks$mean, Models.PQL[[3]]$risks$mean,
     main="PQL. Model ST2 vs. Model ST3", cex.main=1.5, xlab="Model ST2", ylab="Model ST3", cex.lab=1.5)
lines(c(0,3),c(0,3))

plot(Models.PQL[[2]]$risks$mean, Models.PQL[[4]]$risks$mean,
     main="PQL. Model ST2 vs. Model ST4", cex.main=1.5, xlab="Model ST2", ylab="Model ST4", cex.lab=1.5)
lines(c(0,3),c(0,3))

plot(Models.PQL[[3]]$risks$mean, Models.PQL[[4]]$risks$mean,
     main="PQL. Model ST3 vs. Model ST4", cex.main=1.5, xlab="Model ST3", ylab="Model ST4", cex.lab=1.5)
lines(c(0,3),c(0,3))

dev.off()


#########################################################################################################
## Figure 2: Maps of posterior spatial patterns (top row) and posterior temporal patterns (middle row) ##
##           obtained with models ST2, ST3, and ST4. Posterior spatio-temporal patterns (bottom row)   ##
##           obtained with Models ST3 and ST4 are shown for three districts (Agra, Balrampu and Gautam ##
##           Buddha Nagar). Results are from the INLA fit.                                             ##
#########################################################################################################

#################################################
## Spatial patterns (posterior mean estimates) ##
#################################################
xi <- data.frame(ID_area=1:S,
                 M2=unlist(lapply(Models.INLA[[2]]$marginals.lincomb.derived[2:(S+1)], function(x) inla.emarginal(exp,x))),
                 M3=unlist(lapply(Models.INLA[[3]]$marginals.lincomb.derived[2:(S+1)], function(x) inla.emarginal(exp,x))),
                 M4=unlist(lapply(Models.INLA[[4]]$marginals.lincomb.derived[2:(S+1)], function(x) inla.emarginal(exp,x))))

Carto.xi <- merge(carto_up,xi)

paleta <- brewer.pal(9,"YlOrRd")
values <- c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.15)

tmap_mode("plot")
Mapa_xi <- tm_shape(Carto.xi) +
  tm_polygons(col=c("M2","M3","M4"), palette=paleta, title="", legend.show=T,
              legend.reverse=T, style="fixed", breaks=values, interval.closure="left") +
  tm_layout(main.title="", main.title.position="center", legend.text.size=1,
            panel.labels=c("Model ST2","Model ST3","Model ST4"), panel.label.bg.color="lightskyblue",
            legend.outside=T, legend.outside.position="right", legend.frame=F) +
  tm_facets(ncol=3, nrow=1)

print(Mapa_xi)
tmap_save(Mapa_xi, file="figures/Figure2a.pdf", width=12, height=4)


###########################################################
## Spatial patterns (posterior exceedence probabilities) ##
###########################################################
prob.xi <- data.frame(ID_area=1:S,
                      M2=unlist(lapply(Models.INLA[[2]]$marginals.lincomb.derived[2:(S+1)], function(x){1-inla.pmarginal(0, x)})),
                      M3=unlist(lapply(Models.INLA[[3]]$marginals.lincomb.derived[2:(S+1)], function(x){1-inla.pmarginal(0, x)})),
                      M4=unlist(lapply(Models.INLA[[4]]$marginals.lincomb.derived[2:(S+1)], function(x){1-inla.pmarginal(0, x)})))

Carto.prob.xi <- merge(carto_up,prob.xi)

paleta <- brewer.pal(7,"Blues")[-c(3,5)]
values <- c(0,0.1,0.2,0.8,0.9,1)

tmap_mode("plot")
Mapa_prob.xi <- tm_shape(Carto.prob.xi) +
  tm_polygons(col=c("M2","M3","M4"), palette=paleta, title="", legend.show=T,
              legend.reverse=T, style="fixed", breaks=values, interval.closure="left",
              labels=c("[0-0.1)","[0.1-0.2)","[0.2-0.8)","[0.8-0.9)","[0.9-1]")) +
  tm_layout(main.title="", main.title.position="center", legend.text.size=1,
            panel.labels=c("Model ST2","Model ST3","Model ST4"), panel.label.bg.color="lightskyblue",
            legend.outside=T, legend.outside.position="right", legend.frame=F) +
  tm_facets(ncol=3, nrow=1)

print(Mapa_prob.xi)


######################
## Temporal pattern ##
######################
SMR.year <- apply(matrix(Data$O,S,T,byrow=F),2,sum)/apply(matrix(Data$E,S,T,byrow=F),2,sum)

graphics.off()
pdf("figures/Figure2b.pdf", width=12, height=4)

par(mfrow=c(1,3), pty="s")
x <- 1:T

for(i in 2:4){
  temporal <- unlist(lapply(Models.INLA[[i]]$marginals.lincomb.derived[(S+2):(1+S+T)], function(x) inla.emarginal(exp,x)))
  aux <- lapply(Models.INLA[[i]]$marginals.lincomb.derived[(S+2):(1+S+T)], function(x) inla.tmarginal(exp,x))
  q1 <- unlist(lapply(aux, function(x) inla.qmarginal(0.025,x)))
  q2 <- unlist(lapply(aux, function(x) inla.qmarginal(0.975,x)))

  plot(range(x),c(0.7,1.3),type="n",xlab="",ylab="", xaxt="n", main=paste("Model ST",i,sep=""))
  axis(1, at=seq(1,T,4), labels=seq(t.from,t.to,4), las=0)
  X.Vec <- c(x, tail(x, 1), rev(x), x[1])
  Y.Vec <- c(q1, tail(q2, 1), rev(q2), q1[1])
  polygon(X.Vec, Y.Vec, col = "grey", border = NA)
  lines(temporal)
  lines(SMR.year, col="red", lwd=2)
  abline(h=1,lty=2)
}

dev.off()


#############################
## Spatio-temporal pattern ##
#############################
loc <- which(data$dist[1:S] %in% c("Agra","Balrampur","Gautam Buddha Nagar"))

M3.delta <- unlist(lapply(Models.INLA[[3]]$marginals.lincomb.derived[(S+T+2):(1+S+T+S*T)], function(x) inla.emarginal(exp,x)))
M3.aux <- lapply(Models.INLA[[3]]$marginals.lincomb.derived[(S+T+2):(1+S+T+S*T)], function(x) inla.tmarginal(exp,x))
M3.q1 <- unlist(lapply(M3.aux, function(x) inla.qmarginal(0.025,x)))
M3.q2 <- unlist(lapply(M3.aux, function(x) inla.qmarginal(0.975,x)))

M3.Delta <- matrix(M3.delta,S,T,byrow=F)
M3.Q1 <- matrix(M3.q1,S,T,byrow=F)
M3.Q2 <- matrix(M3.q2,S,T,byrow=F)

M4.delta <- unlist(lapply(Models.INLA[[4]]$marginals.lincomb.derived[(S+T+2):(1+S+T+S*T)], function(x) inla.emarginal(exp,x)))
M4.aux <- lapply(Models.INLA[[4]]$marginals.lincomb.derived[(S+T+2):(1+S+T+S*T)], function(x) inla.tmarginal(exp,x))
M4.q1 <- unlist(lapply(M4.aux, function(x) inla.qmarginal(0.025,x)))
M4.q2 <- unlist(lapply(M4.aux, function(x) inla.qmarginal(0.975,x)))

M4.Delta <- matrix(M4.delta,S,T,byrow=F)
M4.Q1 <- matrix(M4.q1,S,T,byrow=F)
M4.Q2 <- matrix(M4.q2,S,T,byrow=F)

color1 <- rgb(red=105,green=105,blue=105,alpha=150,maxColorValue=255)
color2 <- rgb(red=0,green=191,blue=255,alpha=150,maxColorValue=255)

graphics.off()
pdf("figures/Figure2c.pdf", width=12, height=4)

par(mfrow=c(1,3), pty="m")
x <- 1:T

for(i in loc){
  plot(range(x),c(0.4,1.8),type="n",xlab="",ylab="", xaxt="n", main=data$dist[i])
  axis(1, at=seq(1,T,4), labels=seq(t.from,t.to,4), las=0)
  
  X.Vec <- c(x, tail(x, 1), rev(x), x[1])
  Y.Vec <- c(M3.Q1[i,], tail(M3.Q2[i,], 1), rev(M3.Q2[i,]), M3.Q1[i,1])
  polygon(X.Vec, Y.Vec, col=color1, border=NA)
  
  X.Vec <- c(x, tail(x, 1), rev(x), x[1])
  Y.Vec <- c(M4.Q1[i,], tail(M4.Q2[i,], 1), rev(M4.Q2[i,]), M4.Q1[i,1])
  polygon(X.Vec, Y.Vec, col=color2, border=NA)
  
  lines(M3.Delta[i,], lwd=2)
  lines(M4.Delta[i,], col="blue", lwd=2)
  
  abline(h=1, lty=2)
  
  legend("topright", legend=c("Model ST3","Model ST4"), col=c("black","blue"), lwd=2)
}

dev.off()


##################################################################################################
## Figure 3: Final risk estimates obtained with models ST3 and ST4 and INLA in three districts, ##
##           Agra, Balrampu and Gautam Buddha Nagar. Black lines and grey credible intervals    ##
##           corresponds to Model ST3, blue lines and credible intervals to Model ST4.          ##
##           Red lines represent the crude standardized mortality ratios.                       ##
##################################################################################################
loc <- which(data$dist[1:S] %in% c("Agra","Balrampur","Gautam Buddha Nagar"))

SMR <- matrix(Data$O/Data$E,S,T,byrow=F)

M3.risk <- matrix(Models.INLA[[3]]$summary.fitted.values$mean[1:(S*T)], S, T, byrow=F)
M3.q1 <- matrix(Models.INLA[[3]]$summary.fitted.values$`0.025quant`[1:(S*T)], S, T, byrow=F)
M3.q2 <- matrix(Models.INLA[[3]]$summary.fitted.values$`0.975quant`[1:(S*T)], S, T, byrow=F)

M4.risk <- matrix(Models.INLA[[4]]$summary.fitted.values$mean[1:(S*T)], S, T, byrow=F)
M4.q1 <- matrix(Models.INLA[[4]]$summary.fitted.values$`0.025quant`[1:(S*T)], S, T, byrow=F)
M4.q2 <- matrix(Models.INLA[[4]]$summary.fitted.values$`0.975quant`[1:(S*T)], S, T, byrow=F)

color1 <- rgb(red=105,green=105,blue=105,alpha=150,maxColorValue=255)
color2 <- rgb(red=0,green=191,blue=255,alpha=150,maxColorValue=255)

graphics.off()
pdf("figures/Figure3.pdf", width=12, height=4)

par(mfrow=c(1,3), pty="m")
x <- 1:T

for(i in loc){
  plot(range(x),c(0.2,2.8),type="n",xlab="",ylab="", xaxt="n", main=data$dist[i])
  axis(1, at=seq(1,T,4), labels=seq(t.from,t.to,4), las=0)
  
  X.Vec <- c(x, tail(x, 1), rev(x), x[1])
  Y.Vec <- c(M3.q1[i,], tail(M3.q2[i,], 1), rev(M3.q2[i,]), M3.q1[i,1])
  polygon(X.Vec, Y.Vec, col=color1, border=NA)
  
  X.Vec <- c(x, tail(x, 1), rev(x), x[1])
  Y.Vec <- c(M4.q1[i,], tail(M4.q2[i,], 1), rev(M4.q2[i,]), M4.q1[i,1])
  polygon(X.Vec, Y.Vec, col=color2, border=NA)
  
  lines(M3.risk[i,], lwd=2)
  lines(M4.risk[i,], col="blue", lwd=2)
  lines(SMR[i,], col="red", lwd=2)
  
  abline(h=1, lty=2)
  
  legend("topright", legend=c("Model ST3","Model ST4","SMR"), col=c("black","blue","red"), lwd=2)
}

dev.off()


############################################################################################
## Figure B1: Boxplots of correlations between the covariates and the spatial eigenvector ##
##            U_xi69 for each year (left) and correlations between the covariates and the ##
##            temporal eigenvector U_gamma13 for each year (right).                       ##
############################################################################################
eigen.spatial <- eigen(Qs)$vectors[,S-1]
eigen.temporal <- eigen(Qt)$vectors[,T-1]
eigen.ST <- eigen(Qst)$vectors[,S*T-S-T+1]

spatial.eigencor.X <- vector("list",6)
temporal.eigencor.X <- vector("list",6)

for(k in 1:6){
  
  ## Compute (spatial) correlations by year ##
  spatial.cor <- numeric(T)
  
  for(i in 1:T){
    X.year <- Data[Data$ID.year==i, which(names(Data)==paste("X",k,sep=""))]
    spatial.cor[i] <- cor(X.year,eigen.spatial)
  }

  ## Compute (temporal) correlations by region ##
  temporal.cor <- numeric(S)
  
  for(i in 1:S){
    X.area <- Data[Data$ID.area==i, which(names(Data)==paste("X",k,sep=""))]
    temporal.cor[i] <- cor(X.area,eigen.temporal)
  }

  spatial.eigencor.X[[k]] <- spatial.cor
  temporal.eigencor.X[[k]] <- temporal.cor
}


graphics.off()
pdf("figures/FigureB1.pdf", width=12, height=8)

par(mfrow=c(1,2), pty="s")
boxplot(spatial.eigencor.X, main=expression(paste("Correlation between X and ",U[xi[69]])), cex.main=1.5, ylim=c(-1,1), xaxt="n")
axis(1, at=1:k, labels=var.name, las=2)

boxplot(temporal.eigencor.X, main=expression(paste("Correlation between X and ",U[gamma[13]])), cex.main=1.5, ylim=c(-1,1), xaxt="n")
axis(1, at=1:k, labels=var.name, las=2)

dev.off()
