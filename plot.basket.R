#### Plotting der Basket Call Option ####
require(qrng)
require(copula)
require(stringr)
require(graphics)
require(latex2exp)

dvektor <- list(5, 10, 20)
Nvektor <- seq(10000, 200000, by=5000)
### Fuer Clayton Copula ###
## Einlesen der '.csv'-Dateien
## hier nur fuer den Fall Tau = 0.2
## fuer Tau = 0.5 '.tau0.2' in den Dateinamen entfernen
res <- read.csv("res.var.C.CDM..tau0.2.csv")
resvar.C.MO.result <- read.csv("res.var.C.MO..tau0.2.csv")
resvar.C.CDMresult <- read.csv("res.var.C.CDM.tau0.2.csv")
resvar.C.CDM..result <- read.csv("res.var.C.CDM...tau0.2.csv")
resvar.C.MO..result <- read.csv("res.var.C.MO...tau0.2.csv")

for(rowindex in 1:length(dvektor)){
    dimens <- dvektor[[rowindex]]
    # fuer den Fall Tau = 0.5 ersetze '02' in filename durch '05'
    filename <- "Var.Q.Basket.Clayton.CDM.dDIMENSION.Tau02.pdf"
    filename <- str_replace(filename, "DIMENSION", toString(dimens))
    pdf(file=filename)
    # fuer den Fall Tau = 0.5 ersetze '0.2' in label durch '0.5'
    label <- "Varianzsch\344tzer f\374r d = DIMENSION, $\\tau$ = 0.2 "
    label <- str_replace(label, "DIMENSION", toString(dimens))
    plot(Nvektor,res[rowindex,], type="l", log="xy",
    ylim=c(0.000001,150), ylab=TeX(label) ,
    xlab="Stichprobengr\366\337e n",
    main="Clayton Copula, \n Basket Call Option mit Strike = 0.9",
    cex.main="1.5" , mgp = c(2.3,0.75,0) , cex.axis=1 ,
    cex.lab= 1.5, lwd=2, col="chartreuse3") # Halton Folge CDM
    lines(Nvektor,resvar.C.MO.result[rowindex,],
    lty=2, lwd=2) # Halton Folge MO
    lines(Nvektor,resvar.C.CDMresult[rowindex,],
    col="darkorange", lwd=2) # PRN CDM
    lines(Nvektor,resvar.C.CDM..result[rowindex,],
    col="blue", lwd=2) # Sobol Folge CDM
    lines(Nvektor,resvar.C.MO..result[rowindex,],
    col="red", lty=2, lwd=2) # Sobol Folge MO
    lines(Nvektor, 1800/Nvektor, type="b", pch=0 ,
    lty = 1, col="black", lwd=1)
    lines(Nvektor, 900/Nvektor^(1.5), type="b",
    pch=1, lty = 1, col="black" , lwd=1)
    legend( 80000, 100, legend=c("CDM(Ghalton)", "CDM(Sobol)",
    "MO(Ghalton)","MO(Sobol)",  "CDM(PRNG)","1800/n","900/n^(1.5)"),
    lty=c(1, 1, 2, 2 ,1, 1,1,0,0 ),
    pch=c(NA,NA,NA,NA,NA,0,1),
    lwd=c(2, 2, 2, 2, 2, NA,NA), bty ="n",
    col=c("chartreuse3", "blue", "black",
    "red", "darkorange","black","black" ))
}
### Fuer t-Copula ###
## Einlesen der '.csv'-Dateien
## hier nur fuer den Fall Tau = 0.2
## fuer Tau = 0.5 '.tau0.2' in den Dateinamen entfernen
resvar.t.CDMresult <- read.csv("res.var.t.CDM.tau0.2.csv")
resvar.t.CDM.result <- read.csv("res.var.t.CDM..tau0.2.csv")
resvar.t.CDM..result <- read.csv("res.var.t.CDM...tau0.2.csv")

for(rowindex in 1:length(dvektor)){
    dimens <- dvektor[[rowindex]]
    # fuer den Fall Tau = 0.5 ersetze '02' in filename durch '05'
    filename <- "Var.Q.Basket.t.CDM.dDIMENSION.Tau02.pdf"
    filename <- str_replace(filename, "DIMENSION", toString(dimens))
    pdf(file=filename)
    # fuer den Fall Tau = 0.5 ersetze '0.2' in label durch '0.5'
    label <- "Varianzsch\344tzer f\374r d = DIMENSION, $\\tau$ = 0.2 "
    label <- str_replace(label, "DIMENSION", toString(dimens))
    plot(Nvektor,resvar.t.CDMresult[rowindex,], type="l", log="xy",
    ylim=c(0.000001,150), col="darkorange",
    xlab="Stichprobengr\366\337e n",
    ylab=TeX(label),
    main="t-Copula, \n Basket Call Option mit Strike = 0.9",
    cex.main="1.5" , mgp = c(2.3,0.75,0) , cex.axis=1 ,
    cex.lab= 1.5, lty=1, lwd=2)  # PRN
    lines(Nvektor,resvar.t.CDM.result[rowindex,],
    col="chartreuse3", lty=1, lwd=2) # Halton Folge
    lines(Nvektor,resvar.t.CDM..result[rowindex,],
    col="blue", lty=1, lwd=2) # Sobol Folge

    lines(Nvektor, 3000/Nvektor, type="b", pch=0,
    col="black", lty=1, lwd=1)
    lines(Nvektor, 3200/Nvektor^(1.5), type="b",
    pch=1, col="black", lty=1, lwd=1)
    legend( 80000, 100, legend=c("CDM(PRNG)", "CDM(Ghalton)",
    "CDM(Sobol)", "3000/n", "3200/n^1.5"),
    bty ="n", lty=c(1,1,1,0,0),
    pch=c(NA,NA,NA,0,1),
    lwd=c(2, 2, 2, NA,NA),
    col=c("darkorange", "chartreuse3", "blue", "black", "black")
    )
}

