rm(list=ls())

library(INLA)
library(sf)
library(spdep)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))


############################################################################
## Read the dowry death mortality data (or replace it with your own data) ##
############################################################################
data <- read.table(file="../data/DowryDeaths_UttarPradesh.txt", header=TRUE)
str(data)

## Define 'V.area' and 'V.year' variables ##
V.area <- "dist"
V.year <- "year"

S <- length(unique(data[,V.area]))
T <- length(unique(data[,V.year]))

t.from <- min(data[,V.year])
t.to <- max(data[,V.year])


#############################################################
## Prepare the data for INLA models (SPATIO-TEMPORAL CASE) ##
#############################################################

## Order the data by year and area ##
data <- data[order(data[,V.year],data[,V.area]),]
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

ones.S <- matrix(1,S,1)
ones.T <- matrix(1,T,1)
ones.ST <- matrix(1,S*T,1)

p <- 6 # Number of covariates

Beta.df <- as.matrix(Data[,paste("X",1:p,sep="")])


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


###############################
## Read the cartography file ##
###############################
carto_up <- st_read("../data/carto_up")
plot(carto_up$geometry, axes=T)

## Add the neighbourhood graph ##
W <- -Qs
diag(W) <- 0
carto.nb <- mat2listw(W, style="B")$neighbours
plot(carto.nb, st_coordinates(st_centroid(carto_up)), pch=19, cex=0.8, col="blue", add=TRUE)


#########################################
## Define the hyperprior distributions ##
#########################################
sdunif="expression:
  logdens=-log_precision/2;
  return(logdens)"


#############################################
## Compute posterior pattern distributions ##
#############################################
compute.patterns <- TRUE  ## Set compute.patterns=FALSE if posterior patterns are not required

if(compute.patterns){
  source("posterior_lincombs_APredictor.R")
  all.lc.Apredictor <- all.lc
  
  source("posterior_lincombs.R")
  all.lc <- all.lc
}else{
  all.lc <- NULL
  all.lc.Apredictor <- NULL
}


#########################################################################################################
#########################################################################################################
## IMPORTANT NOTE: Uncomment and replace formula lines (f.M1-f.M4) for automatic definition when using ##
##                 your own data. In that case, the covariates should be named as X1, X2, ..., Xp.     ##                  ##
#########################################################################################################
#########################################################################################################
strategy <- "simplified.laplace"

#########################################
## Model ST1: Intercept + fixed-effect ##
#########################################
f.M1 <- O ~ X1 + X2 + X3 + X4 + X5 + X6
# f.M1 <- formula(paste("O ~",paste0(colnames(Data)[grep("^X", colnames(Data))],collapse="+")))

Model1 <- inla(f.M1, family="poisson", data=Data, E=E,
               control.predictor=list(compute=TRUE, cdf=c(log(1))),
               control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
               control.inla=list(strategy=strategy))

Model1$summary.fixed


#########################################################################
## Model ST2: Intercept + random-effects (Eq. 2.4 -> with constraints) ##
#########################################################################

## Sum-to-zero constraints for the interaction ##
R <- kronecker(Qt,Qs)
r.def <- S+T-1
Bst1 <- kronecker(t(ones.T),diag(S))
Bst2 <- kronecker(diag(T),t(ones.S))
Bst <- rbind(Bst1[-1,],Bst2[-1,])

## INLA model ##
f.M2 <- O ~ X1 + X2 + X3 + X4 + X5 + X6 + 
  f(ID.area, model="besag", graph=Qs, constr=TRUE, rankdef=1, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.year, model="rw1", constr=TRUE, rankdef=1, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.area.year, model="generic0", Cmatrix=R, rankdef=r.def, hyper=list(prec=list(prior=sdunif)),
    constr=TRUE, extraconstr=list(A=Bst, e=rep(0,nrow(Bst))))

# f.M2 <- formula(paste("O ~",paste0(colnames(Data)[grep("^X", colnames(Data))],collapse="+"),
#                       "+ f(ID.area, model='besag', graph=Qs, constr=TRUE, rankdef=1, hyper=list(prec=list(prior=sdunif)))",
#                       "+ f(ID.year, model='rw1', constr=TRUE, rankdef=1, hyper=list(prec=list(prior=sdunif)))",
#                       "+ f(ID.area.year, model='generic0', Cmatrix=R, rankdef=r.def, hyper=list(prec=list(prior=sdunif)),
#                            constr=TRUE, extraconstr=list(A=Bst, e=rep(0,nrow(Bst))))"))

Model2 <- inla(f.M2, family="poisson", data=Data, E=E,
               control.predictor=list(compute=TRUE, cdf=c(log(1))),
               control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
               lincomb=all.lc,
               control.inla=list(strategy=strategy))

Model2$summary.fixed


#####################################################################
## Model ST3a: Restricted regression (Eq. 3.2 -> with constraints) ##
#####################################################################
W <- diag(Model2$summary.fitted.values$mode*Data$E)
W.sqrt <- diag(sqrt(diag(W)))

X <- cbind(ones.ST,as.matrix(Beta.df))
P <- W.sqrt%*%X%*%solve(t(X)%*%W%*%X)%*%t(X)%*%W.sqrt
#rankMatrix(P)
Pc <- diag(S*T)-P
#rankMatrix(Pc)

eigen.Pc <- eigen(Pc)
L <- eigen.Pc$vectors[,eigen.Pc$values>1e-12]
# eigen.P <- eigen(P)
# K <- as.matrix(eigen.P$vectors[,eigen.P$values>1e-12])
#rankMatrix(rbind(X,K))[1]==ncol(X)  # K is proportional to X
#all(t(K)%*%L<1e-12)		# K'L=0

M <- solve(W.sqrt)%*%L%*%t(L)%*%W.sqrt
Z.area <- M%*%kronecker(ones.T,diag(S))
Z.year <- M%*%kronecker(diag(T),ones.S)
Z.area.year <- M%*%diag(S*T)

Data.M3a <- list(O=Data$O, E=Data$E,
                 X1=Data$X1, X2=Data$X2, X3=Data$X3, X4=Data$X4, X5=Data$X5, X6=Data$X6,
                 Intercept=c(1,rep(NA,6+S+T+S*T)),
                 beta1=c(NA,1,rep(NA,5+S+T+S*T)),
                 beta2=c(rep(NA,2),1,rep(NA,4+S+T+S*T)),
                 beta3=c(rep(NA,3),1,rep(NA,3+S+T+S*T)),
                 beta4=c(rep(NA,4),1,rep(NA,2+S+T+S*T)),
                 beta5=c(rep(NA,5),1,rep(NA,1+S+T+S*T)),
                 beta6=c(rep(NA,6),1,rep(NA,S+T+S*T)),
                 ID.area=c(rep(NA,1+6),1:S,rep(NA,T+S*T)),
                 ID.year=c(rep(NA,1+6+S),1:T,rep(NA,S*T)),
                 ID.area.year=c(rep(NA,1+6+S+T),1:(S*T)))

# Data.M3a <- list(O=Data$O, E=Data$E)
# Data.M3a <- append(Data.M3a, setNames(lapply(colnames(Beta.df), function(x) Beta.df[,x]), colnames(Beta.df)))
# Data.M3a <- append(Data.M3a, list(Intercept=c(1,rep(NA,p+S+T+S*T))))
# for(i in 1:p){
#    aux <- rep(NA,1+p+S+T+S*T)
#    aux[1+i] <- 1
#    Data.M3a <- append(Data.M3a, setNames(list(aux), paste0("beta",i)))
# }
# Data.M3a <- append(Data.M3a, list(ID.area=c(rep(NA,1+p),1:S,rep(NA,T+S*T))))
# Data.M3a <- append(Data.M3a, list(ID.year=c(rep(NA,1+p+S),1:T,rep(NA,S*T))))
# Data.M3a <- append(Data.M3a, list(ID.area.year=c(rep(NA,1+p+S+T),seq(1,S*T))))


## Sum-to-zero constraints for the interaction ##
R <- kronecker(Qt,Qs)
r.def <- S+T-1
Bst1 <- kronecker(t(ones.T),diag(S))
Bst2 <- kronecker(diag(T),t(ones.S))
Bst <- rbind(Bst1[-1,],Bst2[-1,])


## INLA model ##
f.M3a <- O ~ -1 + Intercept + beta1 +  beta2 + beta3 + beta4 + beta5 + beta6 +
  f(ID.area, model="generic0", Cmatrix=Qs, rankdef=1, constr=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.year, model="generic0", Cmatrix=Qt, rankdef=1, constr=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.area.year, model="generic0", Cmatrix=R, rankdef=r.def, constr=TRUE, hyper=list(prec=list(prior=sdunif)),
    extraconstr=list(A=Bst, e=rep(0,nrow(Bst))))

# f.M3a <- formula(paste("O ~ -1 + Intercept +",paste0('beta',1:p,collapse="+"),
#                        "+ f(ID.area, model='generic0', Cmatrix=Qs, rankdef=1, constr=TRUE, hyper=list(prec=list(prior=sdunif)))",
#                        "+ f(ID.year, model='generic0', Cmatrix=Qt, rankdef=1, constr=TRUE, hyper=list(prec=list(prior=sdunif)))",
#                        "+ f(ID.area.year, model='generic0', Cmatrix=R, rankdef=r.def, hyper=list(prec=list(prior=sdunif)),
#                             constr=TRUE, extraconstr=list(A=Bst, e=rep(0,nrow(Bst))))"))

Apredictor <- cbind(rep(1,S*T),Beta.df,Z.area,Z.year,Z.area.year)

Model3a <- inla(f.M3a, family="poisson", data=Data.M3a, E=E,
                control.predictor=list(compute=TRUE, A=Apredictor),
                control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                control.inla=list(strategy=strategy),
                lincomb=all.lc.Apredictor,
                control.mode=list(theta=Model2$mode$theta, restart=FALSE),
                control.fixed=list(prec=0))

Model3a$summary.fixed


####################################################################
## Model ST3b: Restricted regression (Eq. 3.1 -> reparameterized) ##
####################################################################
W <- diag(Model2$summary.fitted.values$mode*Data$E)
W.sqrt <- diag(sqrt(diag(W)))

X <- cbind(ones.ST,as.matrix(Beta.df))
P <- W.sqrt%*%X%*%solve(t(X)%*%W%*%X)%*%t(X)%*%W.sqrt
#rankMatrix(P)
Pc <- diag(S*T)-P
#rankMatrix(Pc)

eigen.Pc <- eigen(Pc)
L <- eigen.Pc$vectors[,eigen.Pc$values>1e-12]
# eigen.P <- eigen(P)
# K <- as.matrix(eigen.P$vectors[,eigen.P$values>1e-12])
#rankMatrix(rbind(X,K))[1]==ncol(X)  # K is proportional to X
#all(t(K)%*%L<1e-12)		# K'L=0

eigen.Qs <- eigen(Qs)
Us <- eigen.Qs$vectors[,eigen.Qs$values>1e-12]
Ds <- diag(eigen.Qs$values[eigen.Qs$values>1e-12])

eigen.Qt <- eigen(Qt)
Ut <- eigen.Qt$vectors[,eigen.Qt$values>1e-12]
Dt <- diag(eigen.Qt$values[eigen.Qt$values>1e-12])

Ust <- kronecker(Ut,Us)
Dst <- kronecker(Dt,Ds)

M <- solve(W.sqrt)%*%L%*%t(L)%*%W.sqrt
Z.area <- M%*%kronecker(ones.T,Us)
Z.year <- M%*%kronecker(Ut,ones.S)
Z.area.year <- M%*%Ust

Data.M3b <- list(O=Data$O, E=Data$E,
                 X1=Data$X1, X2=Data$X2, X3=Data$X3, X4=Data$X4, X5=Data$X5, X6=Data$X6,
                 Intercept=c(1,rep(NA,6+S-1+T-1+(S-1)*(T-1))),
                 beta1=c(NA,1,rep(NA,5+S-1+T-1+(S-1)*(T-1))),
                 beta2=c(rep(NA,2),1,rep(NA,4+S-1+T-1+(S-1)*(T-1))),
                 beta3=c(rep(NA,3),1,rep(NA,3+S-1+T-1+(S-1)*(T-1))),
                 beta4=c(rep(NA,4),1,rep(NA,2+S-1+T-1+(S-1)*(T-1))),
                 beta5=c(rep(NA,5),1,rep(NA,1+S-1+T-1+(S-1)*(T-1))),
                 beta6=c(rep(NA,6),1,rep(NA,S-1+T-1+(S-1)*(T-1))),
                 ID.area=c(rep(NA,1+6),1:(S-1),rep(NA,T-1+(S-1)*(T-1))),
                 ID.year=c(rep(NA,1+6+S-1),1:(T-1),rep(NA,(S-1)*(T-1))),
                 ID.area.year=c(rep(NA,1+6+S-1+T-1),seq(1,(S-1)*(T-1))))

# Data.M3b <- list(O=Data$O, E=Data$E)
# Data.M3b <- append(Data.M3b, setNames(lapply(colnames(Beta.df), function(x) Beta.df[,x]), colnames(Beta.df)))
# Data.M3b <- append(Data.M3b, list(Intercept=c(1,rep(NA,p+S-1+T-1+(S-1)*(T-1)))))
# for(i in 1:p){
#    aux <- rep(NA,1+p+S-1+T-1+(S-1)*(T-1))
#    aux[1+i] <- 1
#    Data.M3b <- append(Data.M3b, setNames(list(aux), paste0("beta",i)))
# }
# Data.M3b <- append(Data.M3b, list(ID.area=c(rep(NA,1+p),1:(S-1),rep(NA,T-1+(S-1)*(T-1)))))
# Data.M3b <- append(Data.M3b, list(ID.year=c(rep(NA,1+p+S-1),1:(T-1),rep(NA,(S-1)*(T-1)))))
# Data.M3b <- append(Data.M3b, list(ID.area.year=c(rep(NA,1+p+S-1+T-1),seq(1,(S-1)*(T-1)))))


## INLA model ##
f.M3b <- O ~ -1 + Intercept + beta1 +  beta2 + beta3 + beta4 + beta5 + beta6 +
  f(ID.area, model="generic0", Cmatrix=Ds, rankdef=0, constr=FALSE, hyper=list(prec=list(prior=sdunif))) +
  f(ID.year, model="generic0", Cmatrix=Dt, rankdef=0, constr=FALSE, hyper=list(prec=list(prior=sdunif))) +
  f(ID.area.year, model="generic0", Cmatrix=Dst, rankdef=0, constr=FALSE, hyper=list(prec=list(prior=sdunif)))

Apredictor <- cbind(rep(1,S*T),Beta.df,Z.area,Z.year,Z.area.year)

Model3b <- inla(f.M3b, family="poisson", data=Data.M3b, E=E,
                control.predictor=list(compute=TRUE, A=Apredictor),
                control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                control.inla=list(strategy=strategy),
                lincomb=all.lc.Apredictor,
                control.mode=list(theta=Model2$mode$theta, restart=FALSE),
                control.fixed=list(prec=0))

Model3b$summary.fixed


############################################################
## Model ST4: Orthogonality constraints (Eq. 3.3 and 3.4) ##
############################################################
W <- diag(Model2$summary.fitted.values$mode*Data$E)

Bs <- rbind(t(ones.ST)%*%W%*%kronecker(ones.T,diag(S)),
            t(Beta.df)%*%W%*%kronecker(ones.T,diag(S)))

Bt <- rbind(t(ones.ST)%*%W%*%kronecker(diag(T),ones.S),
            t(Beta.df)%*%W%*%kronecker(diag(T),ones.S))

Bst <- rbind(t(Beta.df)%*%W,
             kronecker(t(ones.T),diag(S))%*%W,
             kronecker(diag(T),t(ones.S))%*%W)
Bst <- Bst[-nrow(Bst),]


R <- kronecker(Qt,Qs)
r.def <- S+T-1

f.M4 <- O ~ X1 + X2 + X3 + X4 + X5 + X6 +
  f(ID.area, model="generic0", Cmatrix=Qs, rankdef=7, constr=FALSE, hyper=list(prec=list(prior=sdunif)),
    extraconstr=list(A=Bs, e=rep(0,nrow(Bs)))) +
  f(ID.year, model="generic0", Cmatrix=Qt, rankdef=7, constr=FALSE, hyper=list(prec=list(prior=sdunif)),
    extraconstr=list(A=Bt, e=rep(0,nrow(Bt)))) +
  f(ID.area.year, model="generic0", Cmatrix=R, rankdef=r.def+7, constr=FALSE, hyper=list(prec=list(prior=sdunif)),
    extraconstr=list(A=Bst, e=rep(0,nrow(Bst))))

# f.M4 <- formula(paste("O ~",paste0(colnames(Data)[grep("^X", colnames(Data))],collapse="+"),
#                        "+ f(ID.area, model='generic0', Cmatrix=Qs, rankdef=", p+1, ", constr=FALSE, hyper=list(prec=list(prior=sdunif)), extraconstr=list(A=Bs, e=rep(0,nrow(Bs))))",
#                        "+ f(ID.year, model='generic0', Cmatrix=Qt, rankdef=", p+1, ", constr=FALSE, hyper=list(prec=list(prior=sdunif)), extraconstr=list(A=Bt, e=rep(0,nrow(Bt))))",
#                        "+ f(ID.area.year, model='generic0', Cmatrix=R, rankdef=",r.def+p+1, ", constr=FALSE, hyper=list(prec=list(prior=sdunif)), extraconstr=list(A=Bst, e=rep(0,nrow(Bst))))"))

Model4 <- inla(f.M4, family="poisson", data=Data, E=E,
               control.predictor=list(compute=TRUE, cdf=c(log(1))),
               control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
               lincomb=all.lc,
               control.inla=list(strategy=strategy))

Model4$summary.fixed


######################
## Save the results ##
######################

save(list=c("Model1","Model2","Model3a","Model3b","Model4"), file="DataAnalysis_INLA.Rdata")
