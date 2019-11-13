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

### 1.2. Testfunktion \Psi(\bm(u)) = 3 (u_1^2+...+u_d^2) / d ###
    psi_one <- function(u,dpsi){
        3/dpsi * rowSums(u^2)
    }

### 1.3. Absoluter Fehler ###
abs_err <- function(x){
    abs(mean(x) - 1)
}

### 2. Fall Studie ###
dvektor <- list(5, 15) # Dimensionen
Nvektor <- seq(1e4, 2e5, by=1e4) # Stichprobengroessen
## Initiiere Matrizen zur Speicherung der Resultate ##
ResultErrCCDMp <-matrix(, nrow=length(dvektor), ncol=length(Nvektor))
ResultErrCCDMh <-matrix(, nrow=length(dvektor), ncol=length(Nvektor))
ResultErrCCDMs <-matrix(, nrow=length(dvektor), ncol=length(Nvektor))
ResultErrtCDMp <-matrix(, nrow=length(dvektor), ncol=length(Nvektor))
ResultErrtCDMh <-matrix(, nrow=length(dvektor), ncol=length(Nvektor))
ResultErrtCDMs <-matrix(, nrow=length(dvektor), ncol=length(Nvektor))
ResultErrCMOh <-matrix(, nrow=length(dvektor), ncol=length(Nvektor))
ResultErrCMOs <-matrix(, nrow=length(dvektor), ncol=length(Nvektor))

B <-25 # Anzahl der Wiederholungen

for(rowindex in 1:length(dvektor)){
d <- dvektor[[rowindex]]
for(colindex in 1:length(Nvektor)){
n <- Nvektor[[colindex]]

## Copulas
tau <- 0.2 # Kendall's tau
## Clayton Copula
family.C <- "Clayton"
th.C <- iTau(getAcop(family.C), tau) # Entsprechendes theta
clayton.cop <- onacopulaL(family.C, nacList=list(th.C, 1:d))
## t-Copula mit 3 Freiheitsgraden
family.t <- "t"
nu <- 3 # Freiheitsgrade
th.t <- iTau(ellipCopula(family.t, df=nu), tau) # Entsprechendes rho
t.cop <- ellipCopula(family.t, param=th.t, dim=d, df=nu)


### 2.1. Stichprobengewinnung ###
## Initiiere Vektoren fuer Zwischenresultate
## Fuer Clayton Copula mit CDM und Marshall Olkin
RES.C.CDMp <- vector()
RES.C.CDMh <- vector()
RES.C.CDMs <- vector()
RES.C.MOh <- vector()
RES.C.MOs <- vector()

# Fuer t-Copula mit CDM
RES.t.CDMp <- vector()
RES.t.CDMh <- vector()
RES.t.CDMs <- vector()

for(b in 1:B){
## Gleichverteilte Zufallszahlen fuer CDM
set.seed(270+b)
U.CDMp <- matrix(runif(n*d), ncol=d) # PRN
set.seed(270+b)
U.CDMh <- ghalton(n, d=d) # Halton Folge
set.seed(270+b)
U.CDMs <- sobol(n, d=d, randomize="digital.shift") # Sobol Folge

## Gleichverteilte Zufallszahlen fuer Marshall Olkin
set.seed(270+b)
U.MOh <- ghalton(n, d=d+1) # Halton Folge
set.seed(270+b)
U.MOs <- sobol(n, d=d+1, randomize="digital.shift") # Sobol Folge

## Clayton Copula Stichproben ueber CDM
U.C.CDMp  <- cCopula(U.CDMp,  cop=clayton.cop, inverse=TRUE) # PRN
# Halton Folge
U.C.CDMh <- cCopula(U.CDMh, cop=clayton.cop, inverse=TRUE)
# Sobol Folge
U.C.CDMs <- cCopula(U.CDMs,  cop=clayton.cop, inverse=TRUE)

## Clayton Copula Stichproben ueber Marshall Olkin
# Halton Folge
U.C.MOh <- rClaytonMO(U.MOh, theta=th.C, nMO=n, dMO=d)
# Sobol Folge
U.C.MOs <- rClaytonMO(U.MOs, theta=th.C, nMO=n, dMO=d)

## t-Copula Stichproben ueber CDM
U.t.CDMp  <- cCopula(U.CDMp,  cop=t.cop, inverse=TRUE) # PRN
# Halton Folge
U.t.CDMh <- cCopula(U.CDMh, cop=t.cop, inverse=TRUE)
# Sobol Folge
U.t.CDMs <- cCopula(U.CDMs, cop=t.cop, inverse=TRUE)

### 2.2. Auswertung der Testfunktion ###
## Fuer Clayton Copula
psi.C.CDMp <- psi_one(u=U.C.CDMp, dpsi=d)
psi.C.CDMh <- psi_one(u=U.C.CDMh, dpsi=d)
psi.C.CDMs <- psi_one(u=U.C.CDMs, dpsi=d)
psi.C.MOh <- psi_one(u=U.C.MOh, dpsi=d)
psi.C.MOs <- psi_one(u=U.C.MOs, dpsi=d)

## Fuer t-Copula
psi.t.CDMp <- psi_one(u=U.t.CDMp, dpsi=d)
psi.t.CDMh <- psi_one(u=U.t.CDMh, dpsi=d)
psi.t.CDMs <- psi_one(u=U.t.CDMs, dpsi=d)

## Absoluter Fehler fuer Clayton Copula mit CDM und Marshall Olkin
err.C.CDMp <- abs_err(psi.C.CDMp)
err.C.CDMh <- abs_err(psi.C.CDMh)
err.C.CDMs <- abs_err(psi.C.CDMs)
err.C.MOh <- abs_err(psi.C.MOh)
err.C.MOs <- abs_err(psi.C.MOs)

## Absoluter Fehler fuer t-Copula mit CDM
err.t.CDMp <- abs_err(psi.t.CDMp)
err.t.CDMh <- abs_err(psi.t.CDMh)
err.t.CDMs <- abs_err(psi.t.CDMs)

## Speicherung des absoluten Fehlers fuer jeweils eine Wiederholung "b"
# Fuer Clayton Copula mit CDM und Marshall Olkin
RES.C.CDMp <- c(RES.C.CDMp, err.C.CDMp )
RES.C.CDMh <- c(RES.C.CDMh, err.C.CDMh)
RES.C.CDMs <- c(RES.C.CDMs, err.C.CDMs)
RES.C.MOh <- c(RES.C.MOh, err.C.MOh)
RES.C.MOs <- c(RES.C.MOs, err.C.MOs)

# Fuer t-Copula mit CDM
RES.t.CDMp <- c(RES.t.CDMp , err.t.CDMp )
RES.t.CDMh <- c(RES.t.CDMh , err.t.CDMh )
RES.t.CDMs <- c(RES.t.CDMs , err.t.CDMs )
}

### Speichere den Erwartungswert ueber "B" Wiederholungen
### fuer den absoluten Fehler in die Matrix
### entsprechend der Stichprobengroesse gegeben durch "rowindex"
### und entsprechend der Dimension gegeben durch "colindex"
ResultErrCCDMp[rowindex,colindex]  <- mean(RES.C.CDMp)
ResultErrCCDMh[rowindex,colindex]  <- mean(RES.C.CDMh)
ResultErrCCDMs[rowindex,colindex]  <- mean(RES.C.CDMs)
ResultErrtCDMp[rowindex,colindex]  <- mean(RES.t.CDMp)
ResultErrtCDMh[rowindex,colindex]  <- mean(RES.t.CDMh)
ResultErrtCDMs[rowindex,colindex]  <- mean(RES.t.CDMs)
ResultErrCMOh[rowindex,colindex]  <- mean(RES.C.MOh)
ResultErrCMOs[rowindex,colindex]  <- mean(RES.C.MOs)
}
}
### 3. Speichere Endresultate im ".csv"-Format ab ###
fwrite(ResultErrCCDMp, file="ResultErrCCDMp02.csv")
fwrite(ResultErrCCDMh, file="ResultErrCCDMh02.csv")
fwrite(ResultErrCCDMs, file="ResultErrCCDMs02.csv")
fwrite(ResultErrtCDMp, file="ResultErrtCDMp02.csv")
fwrite(ResultErrtCDMh, file="ResultErrtCDMh02.csv")
fwrite(ResultErrtCDMs, file="ResultErrtCDMs02.csv")
fwrite(ResultErrCMOh, file="ResultErrCMOh02.csv")
fwrite(ResultErrCMOs, file="ResultErrCMOs02.csv")