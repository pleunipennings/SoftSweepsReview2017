#This script will run simulations (using a separate C++ script) to create figure 3B for the Soft Sweeps Review

####################
#Preparations
#setwd("~/Dropbox/SoftSweepsReview/Code/Distributions/")
system("./make_HIV1site_1locus_Distribution")    #compile the C++ code
library(plotrix)

listMutRates<-round((10^-5)/2,7) #only one mutation rate for this figure
options("scipen" = 999)

for (epistasis in c(0)){
    for (sd in c(0.01)){
    for (det in c(0)){
        Ne = 10000 #currently Ne cannot be changed in the sims
        NUMRUNS = 1000;
        sb=0.1
        seed =1
        for (i in 1:length(listMutRates)){
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
    x<-c(x,
	"echo \"", "$seed", "$NUMRUNS", "$mu", "$sd","$sb", "$epi", "$det",
	paste("\" | ./CodeSoftSweepsReview_1locus_Distribution > Dis", det, "u", mu, "_" , sd, "_", "1loc" ,".txt",sep="")
    )
	write(x,file="./tempscript.sh")
	system("chmod 775 ./tempscript.sh")
	#Actually running the simulations takes time. Only needed if you want to recreate the data.
    #system("./tempscript.sh") #Run tempscript.sh 
}}}}

#Read data that were created
    filename<-"Dis0u0.000005_0.01_1loc.txt"
    print(filename)
    x<-read.csv(filename,sep="\t",header = FALSE)

#pdf("NewDistributions.pdf",width=8,height=8) #uncomment if you want to make the figure as a pdf
par(mar = c(4,6,4,1))
hist(x$V4[x$V4>0],xlim=c(0,50),col="darkgrey", ylim=c(0,160),
     breaks=c(seq(0,100,2),1000)-1.5, 
     main = "B. Distribution of number of copies", cex.main=1.8,
     xlab="", ylab="", cex.lab=2,
     yaxt="n", xaxt="n",
     freq=TRUE)
hist(rep(0,160),breaks=c(seq(0,100,2),1000)-1.5,col="darkgrey",ylim=c(0,100),xlim=c(0,50),add=T,freq=TRUE)
axis(2,labels = c(seq(0,100,by=20),660)/1000, at = c(seq(0,100,by=20),160), las=1,cex.axis=1.5)
axis(1,labels = seq(0,50,by=10), at = seq(0,50,by=10),line=-1,cex.axis=1.5)
points(-1:0,c(131,131),pch=15,cex=2.5,col="white")
axis.break(axis=2,breakpos=130,bgcol="white",breakcol="black",style="slash",brw=0.02)
hist(x$V4[x$V2==1],col="#abd9e9",add=T,breaks=c(seq(0,100,2),1000)-1.5,freq=T)
lines(x=c(mu*10000/sd-1.5,mu*10000/sd-1.5),y=c(0,120),lwd=3,lty=2)
mtext("Number of copies",side=1,line=2.3,cex=1.7)
mtext("Frequency", side=2,cex=1.7,line=4)
#dev.off() #uncomment if you want to make the figure as a pdf



