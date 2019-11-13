#### Plotting der Testfunktion ####
require(qrng)
require(copula)
require(stringr)
require(graphics)
require(latex2exp)

dvektor <- list(5, 15) # Dimensionen
Nvektor <- seq(1e4, 2e5, by=1e4) # Stichprobengroessen
### Einlesen der ".csv"-Dateien
ResultErrCCDMp <- read.csv("ResultErrCCDMp02.csv")
ResultErrCCDMh <- read.csv("ResultErrCCDMh02.csv")
ResultErrCCDMs <- read.csv("ResultErrCCDMs02.csv")
ResultErrCMOh <- read.csv("ResultErrCMOh02.csv")
ResultErrCMOs <- read.csv("ResultErrCMOs02.csv")
ResultErrtCDMp <- read.csv("ResultErrtCDMp02.csv")
ResultErrtCDMh <- read.csv("ResultErrtCDMh02.csv")
ResultErrtCDMs <- read.csv("ResultErrtCDMs02.csv")

### Fuer t-Copula ###
for(rowindex in 1:length(dvektor)){
    dimens <- dvektor[[rowindex]]
    filename <- "Err.t.CDM.dDIMENSION.Tau02.pdf"
    filename <- str_replace(filename, "DIMENSION", toString(dimens))
    pdf(file=filename)
    label <- "Absolute Fehler f\374r d = DIMENSION, $\\tau$ = 0.2 "
    label <- str_replace(label, "DIMENSION", toString(dimens))
    plot(Nvektor,ResultErrtCDMp[rowindex,], type="l", log="xy",
    ylim=c(0.000001,150), col="darkorange",
    xlab="Stichprobengr\366\337e n",
    ylab=TeX(label), main="t-Copula", cex.main="1.5" ,
    mgp = c(2.3,0.75,0) , cex.axis=1 ,
    cex.lab= 1.5, lty=1, lwd=2)  # PRN
    lines(Nvektor,ResultErrtCDMh[rowindex,],
    col="chartreuse3", lty=1, lwd=2) # Halton Folge
    lines(Nvektor,ResultErrtCDMs[rowindex,],
    col="blue", lty=1, lwd=2) # Sobol Folge

    lines(Nvektor, 1/Nvektor, type="b", pch=0,
    col="black", lty=1, lwd=1)
    lines(Nvektor, 1/Nvektor^(0.5), type="b",pch=1,
    col="black", lty=1, lwd=1)
    legend( 80000, 100,
    legend=c("CDM(PRNG)", "CDM(Ghalton)",
    "CDM(Sobol)", "1/n", "1/n^0.5"),
    bty ="n", lty=c(1,1,1,0,0),
    pch=c(NA,NA,NA,0,1), lwd=c(2, 2, 2, NA,NA),
    col=c("darkorange", "chartreuse3",
    "blue", "black", "black")
    )
}

### Fuer Clayton Copula ###
for(rowindex in 1:length(dvektor)){
    dimens <- dvektor[[rowindex]]
    filename <- "Err.Clayton.CDM.dDIMENSION.Tau02.pdf"
    filename <- str_replace(filename, "DIMENSION", toString(dimens))
    pdf(file=filename)
    label <- "Absoluter Fehler f\374r d = DIMENSION, $\\tau$ = 0.2 "
    label <- str_replace(label, "DIMENSION", toString(dimens))
    plot(Nvektor,ResultErrCCDMp[rowindex,], type="l", log="xy",
    ylim=c(0.000001,150),
    ylab=TeX(label) ,
    xlab="Stichprobengr\366\337e n",
    main="Clayton Copula",
    cex.main="1.5" , mgp = c(2.3,0.75,0) ,
    cex.axis=1 , cex.lab= 1.5, lwd=2,
    col="darkorange") # PRN CDM
    lines(Nvektor,ResultErrCCDMh[rowindex,],
    col="chartreuse3", lwd=2) # Halton Folge CDM
    lines(Nvektor,ResultErrCCDMs[rowindex,],
    col="blue", lwd=2) # Sobol Folge CDM
    lines(Nvektor,ResultErrCMOh[rowindex,],
    lty=2, lwd=2) # Halton Folge MO
    lines(Nvektor,ResultErrCMOs[rowindex,],
    col="red", lty=2, lwd=2) # Sobol Folge MO
    lines(Nvektor, 1/Nvektor, type="b", pch=0 ,
    lty = 1, col="black", lwd=1)
    lines(Nvektor, 1/Nvektor^(0.5), type="b", pch=1,
    lty = 1, col="black" , lwd=1)
    legend( 80000, 100, legend=c( "CDM(PRNG)", "CDM(Ghalton)",
    "CDM(Sobol)", "MO(Ghalton)","MO(Sobol)","1/n","1/n^(0.5)"),
    lty=c(1, 1,1,2,2,0,0 ), pch=c(NA,NA,NA,NA,NA,0,1),
    lwd=c(2, 2, 2, 2, 2, NA,NA), bty ="n",
    col=c("darkorange", "chartreuse3", "blue",
    "black", "red", "black","black" ))
}