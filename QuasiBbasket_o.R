#### Basierend auf den Demos "basket_options" und ####
#### "test_functions" aus dem qrng Paket ####
require(qrng)
require(copula)
require(stringr)
require(data.table)

### 1. Definiere Funktionen ###
### 1.1. Marshall-Olkin (MO) fuer Clayton Copulas ###
rClaytonMO <- function(u, theta, nMO, dMO)
{
    if(!is.matrix(u)) u <- rbind(u)

    dim. <- dim(u)
    dMO <- dim.[2] - 1
    V <- qgamma(u[,dMO+1], shape=1/theta) # n-Vektor
    E <- qexp(u[,seq_len(dMO)]) # (n,d)-Matrix

    (1 + E/matrix(rep(V, dMO), ncol=dMO))^(-1/theta) # (n,d)-Matrix
}
## Berechnung der Stichproben einer d-dimensionalen
## geometrischen brownsche Bewegung (GBM)
## uebernommen aus "basket_options"
## u (n,d)-Matrix von Stichproben in [0,1]
## S0 d-Vektor
## mu d-Vektor enthaelt den Drift von GBMs
## sigma d-Vektor enthaelt die Volatilitaet von GBMs
rGeoBM <- function(u, S0, mu, sigma, T, dBM, nBM)
{
    stopifnot(0 < u, u < 1, length(mu) == (dBM <- length(S0)), mu >= 0,
              length(sigma) == dBM, sigma >= 0, T > 0)
    log.diff <- qnorm(u) * matrix(rep(sigma, each=nBM), ncol=dBM)
    log.drft <- (mu - sigma^2 / 2) * T # d-Vektor
    log. <- matrix(rep(log.drft, nBM),
    ncol=dBM, byrow=TRUE) + log.diff # (n,d)-Matrix
    matrix(rep(S0, nBM), ncol=dBM,
    byrow=TRUE) * exp(log.) # S_t, t in 1,..,T; (n,d)-Matrix
}
## Auszahlungsfunktion uebernommen aus "basket_options"
## K Strike einer Option
## S0 d-Vektor
## S (n, d)-Matrix
payoff <- function(K, N, S0, S, type = c("call", "put"),
                   method = c("basket", "worst.of", "best.of"))
{
    stopifnot(K >= 0, N >= 0, S0 >= 0, S >= 0, length(S0) == ncol(S))
    type <- match.arg(type)
    method <- match.arg(method)
    perf <- switch(method,
                   "basket" = {
        rowMeans(t(t(S)/S0))
    },
    "worst.of" = {
        apply(t(t(S)/S0), 1, min)
    },
    "best.of" = {
        apply(t(t(S)/S0), 1, max)
    },
    stop("Wrong 'method'"))
    N * pmax(0, if(type=="call") perf - K else K - perf)
}
### 2. Fall Studie ###
dvektor <- list(5, 10, 20) # Dimesionen
Nvektor <- seq(10000, 200000, by=5000) # Stichprobengroesse
## Initiiere Matrizen zur Speicherung der Resultate ##
resvar.C.CDM.result <-matrix(, nrow=length(dvektor),
ncol=length(Nvektor))
resvar.C.CDMresult <-matrix(, nrow=length(dvektor),
ncol=length(Nvektor))
resvar.C.MO.result <-matrix(, nrow=length(dvektor),
ncol=length(Nvektor))
resvar.C.CDM..result <-matrix(, nrow=length(dvektor),
ncol=length(Nvektor))
resvar.C.MO..result <-matrix(, nrow=length(dvektor),
ncol=length(Nvektor))
resvar.t.CDMresult <-matrix(, nrow=length(dvektor),
ncol=length(Nvektor))
resvar.t.CDM.result<-matrix(, nrow=length(dvektor),
ncol=length(Nvektor))
resvar.t.CDM..result <-matrix(, nrow=length(dvektor),
ncol=length(Nvektor))

for(rowindex in 1:length(dvektor)){
d <- dvektor[[rowindex]]
for(colindex in 1:length(Nvektor)){
n <- Nvektor[[colindex]]

### 2.1. Stochastische Parameter ###
sigma <- rep(0.2, d) # Volatilitaet
r <- 0.0001 # Stetig verzinste Short rate
S0 <- rep(100, d) # Anfangsbestand
K <- 1.1 # Strike
N <- 1000
T <- 1 # Zeit

## Copulas
## hier nur fuer den Fall Tau = 0.5
## fuer den Fall Tau = 0.2 setze tau <- 0.2
tau <- 0.5 # Kendall's tau
## Clayton Copula
family.C <- "Clayton"
th.C <- iTau(getAcop(family.C), tau) # Entsprechendes theta
clayton.cop <- onacopulaL(family.C, nacList=list(th.C, 1:d))
## t-Copula mit 3 Freiheitsgraden
family.t <- "t"
nu <- 3 # Freiheitsgrade
th.t <- iTau(ellipCopula(family.t, df=nu), tau) # Entsprechendes rho
t.cop <- ellipCopula(family.t, param=th.t, dim=d, df=nu)
B <-25 # Anzahl der Wiederholungen

### 2.2. Stichprobengewinnung ###
## Initiiere Vektoren fuer Zwischenresultate
## Furr t-Copula mit CDM
RES.t.CDM <- vector()
RES.t.CDM. <- vector()
RES.t.CDM.. <- vector()

## Fuer Clayton Copula  mit Marshall Olkin
RES.C.MO. <- vector()
RES.C.MO.. <- vector()

## Fuer Clayton Copula mit CDM
RES.C.CDM <- vector()
RES.C.CDM. <- vector()
RES.C.CDM.. <- vector()


for(b in 1:B){
## Gleichverteilte Zufallszahlen fuer CDM
set.seed(270+b)
U.CDM  <- matrix(runif(n*d), ncol=d) # PRN
set.seed(270+b)
U.CDM. <- ghalton(n, d=d) # Halton Folge
set.seed(270+b)
U.CDM.. <- sobol(n, d=d, randomize="digital.shift") # Sobol Folge

## Gleichverteilte Zufallszahlen fuer MO
set.seed(270+b)
U.MO. <- ghalton(n, d=d+1) # Halton Folge
set.seed(270+b)
U.MO.. <- sobol(n, d=d+1, randomize="digital.shift") # Sobol Folge

## t-Copula Stichproben ueber CDM
U.t.CDM  <- cCopula(U.CDM,  cop=t.cop, inverse=TRUE) # PRN
# Halton Folge
U.t.CDM. <- cCopula(U.CDM., cop=t.cop, inverse=TRUE)
# Sobol Folge
U.t.CDM.. <- cCopula(U.CDM.., cop=t.cop, inverse=TRUE)

## Clayton Copula Stichproben ueber CDM
U.C.CDM  <- cCopula(U.CDM,  cop=clayton.cop, inverse=TRUE) # PRN
# Halton Folge
U.C.CDM. <- cCopula(U.CDM., cop=clayton.cop, inverse=TRUE)
# Sobol Folge
U.C.CDM.. <- cCopula(U.CDM..,  cop=clayton.cop, inverse=TRUE)

## Clayton Copula Stichproben ueber MO
# Halton Folge
U.C.MO. <- rClaytonMO(U.MO., theta=th.C, nMO=n, dMO=d)
# Sobol Folge
U.C.MO.. <- rClaytonMO(U.MO.., theta=th.C, nMO=n, dMO=d)

## Geometrische brownsche Bewegungen Stichproben
S.t.CDM <- rGeoBM(U.t.CDM, S0=S0, mu=rep(r, d), sigma=sigma,
T=T, nBM=n, dBM=d)
S.C.CDM <- rGeoBM(U.C.CDM, S0=S0, mu=rep(r, d), sigma=sigma,
T=T, nBM=n, dBM=d)

## Quasi-Halton-Geometrische brownsche Bewegungen Stichproben
S.t.CDM. <- rGeoBM(U.t.CDM., S0=S0, mu=rep(r, d), sigma=sigma,
T=T, nBM=n, dBM=d)
S.C.CDM. <- rGeoBM(U.C.CDM., S0=S0, mu=rep(r, d), sigma=sigma,
T=T, nBM=n, dBM=d)
S.C.MO.  <- rGeoBM(U.C.MO.,  S0=S0, mu=rep(r, d), sigma=sigma,
T=T, nBM=n, dBM=d)

## Quasi-Sobol-Geometrische brownsche Bewegungen Stichproben
S.t.CDM.. <- rGeoBM(U.t.CDM.., S0=S0, mu=rep(r, d),
sigma=sigma, T=T, nBM=n, dBM=d)
S.C.CDM.. <- rGeoBM(U.C.CDM.., S0=S0, mu=rep(r, d),
sigma=sigma, T=T, nBM=n, dBM=d)
S.C.MO..  <- rGeoBM(U.C.MO..,  S0=S0, mu=rep(r, d),
sigma=sigma, T=T, nBM=n, dBM=d)

### 2.3. Funktionale Berechnung ###

erT <- exp(-r*T)

## Nutzung der PRN Stichproben
## Basket Call
basket.t.CDM <- erT * mean(payoff(K, N=N, S0=S0, S=S.t.CDM))
basket.C.CDM <- erT * mean(payoff(K, N=N, S0=S0, S=S.C.CDM))

## Nutzung der Halton Stichproben
## Basket Call
basket.t.CDM. <- erT * mean(payoff(K, N=N, S0=S0, S=S.t.CDM.))
basket.C.CDM. <- erT * mean(payoff(K, N=N, S0=S0, S=S.C.CDM.))
basket.C.MO.  <- erT * mean(payoff(K, N=N, S0=S0, S=S.C.MO.))

## Nutzung der Sobol Stichproben
## Basket Call
basket.t.CDM.. <- erT * mean(payoff(K, N=N, S0=S0, S=S.t.CDM..))
basket.C.CDM.. <- erT * mean(payoff(K, N=N, S0=S0, S=S.C.CDM..))
basket.C.MO..  <- erT * mean(payoff(K, N=N, S0=S0, S=S.C.MO..))

#Fuer t-Copula
RES.t.CDM.. <- c(RES.t.CDM.. , basket.t.CDM..) # Sobol Folge
RES.t.CDM. <- c(RES.t.CDM. , basket.t.CDM.) # Halton Folge
RES.t.CDM <- c(RES.t.CDM , basket.t.CDM) # PRN

#Fuer Clayton Copula
RES.C.CDM <- c(RES.C.CDM, basket.C.CDM)
RES.C.CDM. <- c(RES.C.CDM.,basket.C.CDM.)
RES.C.CDM.. <- c(RES.C.CDM..,basket.C.CDM..)
RES.C.MO. <- c(RES.C.MO., basket.C.MO.)
RES.C.MO.. <- c(RES.C.MO.., basket.C.MO..)
}

### Varianzberechnung fuer Clayton Copula mit CDM
var.C.CDM <- var(RES.C.CDM)
var.C.CDM. <- var(RES.C.CDM.)
var.C.CDM.. <- var(RES.C.CDM..)
### Varianzberechnung fuer Clayton Copula mit Marshall
var.C.MO. <- var(RES.C.MO.)
var.C.MO.. <- var(RES.C.MO..)
### Varianzberechnung fuer t-Copula mit CDM
var.t.CDM <- var(RES.t.CDM)
var.t.CDM. <- var(RES.t.CDM.)
var.t.CDM.. <- var(RES.t.CDM..)

### Speichere den Wert fuer die Varianz in die Matrix
### entsprechend der Stichprobengroesse gegeben durch "rowindex"
### und entsprechend der Dimension gegeben durch "colindex"
## fuer Clayton Copula mit CDM
resvar.C.CDMresult[rowindex, colindex] <- var.C.CDM
resvar.C.CDM.result[rowindex,colindex] <- var.C.CDM.
resvar.C.CDM..result[rowindex,colindex] <- var.C.CDM..
## fuer Clayton Copula mit Marshall Olkin
resvar.C.MO.result[rowindex, colindex] <- var.C.MO.
resvar.C.MO..result[rowindex, colindex] <- var.C.MO..
## fuer t Copula mit CDM
resvar.t.CDMresult[rowindex,colindex] <- var.t.CDM
resvar.t.CDM.result[rowindex,colindex] <- var.t.CDM.
resvar.t.CDM..result[rowindex,colindex] <- var.t.CDM..
}
}

### 3. Speichere Endresultate im ".csv"-Format ab ###
## hier nur fuer den Fall Tau = 0.5
## fuer den Fall Tau = 0.2 ersetze '.csv' durch 'tau0.2.csv'
## fuer Clayton Copula mit CDM
fwrite(resvar.C.CDMresult, file="res.var.C.CDM.csv")
fwrite(resvar.C.CDM.result, file="res.var.C.CDM..csv")
fwrite(resvar.C.CDM..result, file="res.var.C.CDM...csv")
## fuer Clayton Copula mit Marshal Olkin
fwrite(resvar.C.MO.result, file="res.var.C.MO..csv")
fwrite(resvar.C.MO..result, file="res.var.C.MO...csv")
## fuer t-Copula mit CDM
fwrite(resvar.t.CDMresult, file="res.var.t.CDM.csv")
fwrite(resvar.t.CDM.result, file="res.var.t.CDM..csv")
fwrite(resvar.t.CDM..result, file="res.var.t.CDM...csv")
