#This script will run simulations (using a separate C++ script) to create figure 3A for the Soft Sweeps Review

####################
#Preparations
#setwd("~/Dropbox/SoftSweepsReview/Code/Det_Vs_Stoch/")
system("./make_HIV1site_1locusDet")    #compile the C++ code

listMutRates<-round(10^(seq(-7,-3.5,by=.5))/2,8)
#listMutRates<-listMutRates[1]
options("scipen" = 999)

epistasis = 0
for (sd in c(0.01)){
    for (det in c(0)){
        NUMRUNS = 2000;
        newmut = 0
        sb=0.1
        seed =1
        for (i in 1:length(listMutRates)){
            if (i==1)NUMRUNS=500000
            if (i==2)NUMRUNS=40000
            if (i==3)NUMRUNS=20000
            if (i==4)NUMRUNS=20000
            if (i>4)NUMRUNS=4000
            mu = listMutRates[i]
            print(paste("mu",mu))
            #make script
            x<-"#!/bin/bash"
            x<-c(x,paste("seed=",seed,sep=""))
            x<-c(x,paste("NUMRUNS=",NUMRUNS,sep=""))
            x<-c(x,paste("sd=",sd,sep=""))
            x<-c(x,paste("sb=",sb,sep=""))    
            x<-c(x,paste("epi=",epistasis,sep=""))
            x<-c(x,paste("mu=",mu,sep=""))
            x<-c(x,paste("det=",det,sep=""))
            x<-c(x,paste("newmut=",newmut,sep=""))
            x<-c(x,
                 "echo \"", "$seed", "$NUMRUNS", "$mu", "$sd","$sb", "$epi", "$det", "$newmut",
                 paste("\" | ./CodeSoftSweepsReview_1locusDet > Det_out_", det, "u", mu, "_" , sd, "_", "1loc" ,".txt",sep="")
            )
            write(x,file="./tempscript.sh")
            system("chmod 775 ./tempscript.sh")
            #Actually running the simulations takes time. Only needed if you want to recreate the data.
            #system("./tempscript.sh") #Run tempscript.sh 
        }
        
        #Read the output of the sims into a dataframe called Data
        Data<-data.frame("N"=0, "Theta_trait"=0,"Theta_locus"=0,"NumHard" = 0,"NumHardSGV" = 0,"NumSos" = 0, "NumMos" = 0, "NoSweep" = 0)
        for (i in 1:length(listMutRates)){
            #for (i in 1:2){
            mu = listMutRates[i]
            print(paste("mu",mu, "epistasis", epistasis))
            filename <- paste("Det_out_",det,"u", mu, "_" , sd, "_", "1loc" ,".txt",sep="")
            print(filename)
            if (file.exists(filename)){
                x<-read.csv(filename,sep="\t",header = FALSE)
                Data[i,]<-x[1,c(2,4,6,8,10,12,14, 16)]
            }
        }    
        Data$FracSoft=(Data$NumSos+Data$NumMos)/(Data$NumSos+Data$NumMos+Data$NumHardSGV)
    }
}

PmultDet  <- function(theta,sd, sb){
    A = 1-(1+theta*sb/sd)*exp(-theta*sb/sd)
    B=1-exp(-theta*sb/sd)
    P = A/B
    return(P)
}

PmultStoch  <- function(theta,sd, sb, N){
    R = sb/(sd+1/(4*N))
    A = theta*R/(1+R)
    B = (1+R)^theta -1
    P = 1- A/B
    return(P)
}

sb=0.1
sd=0.01
thetalist<-10^seq(-3,1,by=0.1)
N=10000

#pdf("DetVsStoch_new.pdf",width=8,height=8) #uncomment if you want to make the figure as a pdf
par(mar = c(4,5,4,1))
plot(log10(thetalist),PmultDet(thetalist,sd,sb),col=1,lty = 2,type="l",lwd=3,
     xlab = "", cex.lab=2,
     ylab= "Probability soft sweep", 
     main="A. Probability of soft sweep from SGV", cex.main=1.8,xaxt="n",yaxt="n")
points(log10(thetalist),PmultStoch(thetalist,sd,sb, N),type="l",lty=1,col=1,lwd=3)
points(log10(Data$Theta_trait), Data$FracSoft,pch=16,cex=1.3)
arrows(log10(Data$Theta_trait), 
       Data$FracSoft-sqrt(Data$FracSoft*(1-Data$FracSoft)/(Data$NumSos+Data$NumMos+Data$NumHardSGV)),
       log10(Data$Theta_trait),
       Data$FracSoft+sqrt(Data$FracSoft*(1-Data$FracSoft)/(Data$NumSos+Data$NumMos+Data$NumHardSGV)),
       code=3, length=0.06, angle = 90,lwd=1.3)
axis(1,10^seq(-3,1,by=0.5), at = log10(10^seq(-3,1,by=0.5)),cex.axis=1.5)
mtext(expression(Theta),side=1,line=2.3,cex=1.7)
axis(2,labels = seq(0,1,by=0.2), at = seq(0,1,by=0.2), las=1,cex.axis=1.5)
abline(v=log10(0.114),lty=3,lwd=3.5)
#dev.off() #uncomment if you want to make the figure as a pdf

