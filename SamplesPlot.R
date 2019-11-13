#### Basierend auf den Demos "basket_options" und ####
#### "test_functions" aus dem qrng Paket ####
require(qrng)
require(copula)
require(stringr)
require(data.table)
require(latex2exp)

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

### 2. Fall Studie ###
d <- 2 # Dimension
n <- 4000 # Stichprobengroesse
N <-1000 # Stichprobengroesse fuer Beispiele
# tcopPRNG1000, tcopQRNG1000, ccopPRNG1000, ccopQRNG1000

## Copulas
tau <- 0.5 # Kendall's tau
## Clayton Copula
family.C <- "Clayton"
th.C <- iTau(getAcop(family.C), tau) # Entsprechendes theta
clayton.cop <- onacopulaL(family.C, nacList=list(th.C, 1:d))
## t-Copula mit 3 Freiheitgraden
family.t <- "t"
nu <- 3 # Freiheitgrade
th.t <- iTau(ellipCopula(family.t, df=nu), tau) # Entsprechendes rho
t.cop <- ellipCopula(family.t, param=th.t, dim=d, df=nu)

### 2.1. Gleichverteilte Stichproben fuer Beispiele
### Quasizufallszahlen vs. Pseudozufallszahlen
## generalisierte Halton Folge
set.seed(271)
Ugh <-  ghalton(N, d=d, method = c("generalized"))
## Sobol Folge ohne Randomisierung
set.seed(271)
Us <-  sobol(N, d=d, randomize = c("none"))
## Sobol Folge mit Randomisierung
set.seed(271)
Urs <- sobol(N, d=d, randomize="digital.shift")
## gleichverteilte Zufallszahlen
set.seed(271)
Upseudo <- matrix(runif(N*d), ncol=d)

### 2.2. Gleichverteilte Stichproben fuer Beispiele
### tcopPRNG1000, tcopQRNG1000, ccopPRNG1000, ccopQRNG1000
set.seed(271)
Up <- matrix(runif(N*d), ncol=d)
set.seed(271)
Uq <-  ghalton(N, d=d)
# t-copula mit PRNG und N=1000
Up.t <- cCopula(Up,  cop=t.cop, inverse=TRUE)
# t-copula mit QRNG und N=1000
Uq.t <- cCopula(Uq,  cop=t.cop, inverse=TRUE)
# clayton copula mit PRNG und N=1000
Up.c <- cCopula(Up,  cop=clayton.cop, inverse=TRUE)
# clayton copula mit QRNG und N=1000
Uq.c <- cCopula(Uq,  cop=clayton.cop, inverse=TRUE)

### 2.3. Gleichverteile Stichproben fuer Beispiele CDM vs. MO
### fuer CDM mit n Stichprobengroesse
set.seed(271)
U.CDM  <- matrix(runif(n*d), ncol=d) # PRN
set.seed(271)
U.CDM. <- ghalton(n, d=d,) # Halton Folge
set.seed(271)
U.CDM.. <- sobol(n, d=d, randomize="digital.shift") # Sobol Folge

### Gleichverteile Stichproben fuer Beispiele CDM vs. MO
### fuer MO mit n Stichprobengroesse
set.seed(271)
U.MO  <- matrix(runif(n*(d+1)), ncol=d+1) # PRN
set.seed(271)
U.MO. <- ghalton(n, d=d+1) # Halton Folge
set.seed(271)
U.MO.. <- sobol(n, d=d+1, randomize=TRUE) # Sobol Folge

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
U.C.MO  <- rClaytonMO(U.MO,  theta=th.C, nMO=n, dMO=d) # PRN
# Halton Folge
U.C.MO. <- rClaytonMO(U.MO., theta=th.C, nMO=n, dMO=d)
# Sobol Folge
U.C.MO.. <- rClaytonMO(U.MO.., theta=th.C, nMO=n, dMO=d)

### 3 PLOT
### 3.1. PLOT fuer Beispiele tcopPRNG1000, tcopQRNG1000,
### ccopPRNG1000, ccopQRNG1000
### und Speicherung der Plots in '.pdf'-Dateien
pdf(file = "tcopPRNG1000.pdf")
plot(Up.t ,xlab=TeX("U_1"), ylab = TeX("U_2") ,
mgp = c(2.3,0.75,0) , cex.axis=2 , cex.lab= 2)
pdf(file = "tcopQRNG1000.pdf")
plot(Uq.t ,xlab=TeX("U_1"), ylab = TeX("U_2") ,
mgp = c(2.3,0.75,0) , cex.axis=2 , cex.lab= 2)
pdf(file = "ccopPRNG1000.pdf")
plot(Up.c ,xlab=TeX("U_1"), ylab = TeX("U_2") ,
mgp = c(2.3,0.75,0) , cex.axis=2 , cex.lab= 2)
pdf(file = "ccopQRNG1000.pdf")
plot(Uq.c ,xlab=TeX("U_1"), ylab = TeX("U_2") ,
mgp = c(2.3,0.75,0) , cex.axis=2 , cex.lab= 2)

### 3.2. PLOT fuer Beispiele CDM vs. MO mit n Stichprobengroesse
### und Speicherung der Plots in '.pdf'-Dateien
pdf(file="Samples.U.t.CDM.pdf")
plot(U.t.CDM, xlab=TeX("U_1"), ylab = TeX("U_2") ,
mgp = c(2.3,0.75,0) , cex.axis=2 , cex.lab= 2,
main="t-Copula; Pseudozufallszahlen; CDM")
pdf(file="Samples.U.t.CDM..pdf")
plot(U.t.CDM., xlab=TeX("U_1"), ylab = TeX("U_2") ,
mgp = c(2.3,0.75,0) , cex.axis=2 , cex.lab= 2,
main="t-Copula; Halton-Folge; CDM")
pdf(file="Samples.U.t.CDM...pdf")
plot(U.t.CDM.., xlab=TeX("U_1"), ylab = TeX("U_2") ,
mgp = c(2.3,0.75,0) , cex.axis=2 , cex.lab= 2,
main="t-Copula; Sobol-Folge; CDM")

pdf(file="Samples.U.C.CDM.pdf")
plot(U.C.CDM, xlab=TeX("U_1"), ylab = TeX("U_2") ,
mgp = c(2.3,0.75,0) , cex.axis=2 , cex.lab= 2,
main="Clayton Copula; Pseudozufallszahlen; CDM")
pdf(file="Samples.U.C.CDM..pdf")
plot(U.C.CDM., xlab=TeX("U_1"), ylab = TeX("U_2") ,
mgp = c(2.3,0.75,0) , cex.axis=2 , cex.lab= 2,
main="Clayton Copula; Halton-Folge; CDM")
pdf(file="Samples.U.C.CDM...pdf")
plot(U.C.CDM.., xlab=TeX("U_1"), ylab = TeX("U_2") ,
mgp = c(2.3,0.75,0) , cex.axis=2 , cex.lab= 2,
main="Clayton Copula; Sobol-Folge; CDM")

pdf(file = "Samples.U.C.MO.pdf")
plot(U.C.MO, xlab=TeX("U_1"), ylab = TeX("U_2") ,
mgp = c(2.3,0.75,0) , cex.axis=2 , cex.lab= 2,
main="Clayton Copula; Pseudozufallszahlen; Marshall Olkin")
pdf(file = "Samples.U.C.MO..pdf")
plot(U.C.MO. , xlab=TeX("U_1"), ylab = TeX("U_2") ,
mgp = c(2.3,0.75,0) , cex.axis=2 , cex.lab= 2,
main="Clayton Copula; Halton-Folge; Marshall Olkin")
pdf(file = "Samples.U.C.MO...pdf")
plot(U.C.MO.. ,xlab=TeX("U_1"), ylab = TeX("U_2") ,
mgp = c(2.3,0.75,0) , cex.axis=2 , cex.lab= 2,
main="Clayton Copula; Sobol-Folge; Marshall Olkin")

### 3.3. PLOT fuer Beispiele
### Quasizufallszahlen vs. Pseudozufallszahlen,
### mit n Stichprobengroesse
### und Speicherung der Plots in '.pdf'-Dateien
pdf(file = "UCDMp.pdf")
plot(Upseudo ,xlab=TeX("U_1"), ylab = TeX("U_2") ,
mgp = c(2.3,0.75,0) , cex.axis=2 , cex.lab= 2,
main="Pseudozufallszahlen")
pdf(file = "UCDMrs.pdf")
plot(Urs ,xlab=TeX("U_1"), ylab = TeX("U_2") ,
mgp = c(2.3,0.75,0) , cex.axis=2 , cex.lab= 2,
main="Sobol-Folge mit digitalem Shift")
pdf(file = "UCDMs.pdf")
plot(Us ,xlab=TeX("U_1"), ylab = TeX("U_2") ,
mgp = c(2.3,0.75,0) , cex.axis=2 , cex.lab= 2,
main="Sobol-Folge")
pdf(file = "UCDMrh.pdf")
plot(Ugh ,xlab=TeX("U_1"), ylab = TeX("U_2") ,
mgp = c(2.3,0.75,0) , cex.axis=2 , cex.lab= 2,
main="generalisierte Halton-Folge")
