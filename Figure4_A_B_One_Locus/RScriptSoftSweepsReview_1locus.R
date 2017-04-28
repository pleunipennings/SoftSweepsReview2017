#This script will run simulations (using a separate C++ script) to create figure 4A, B for the Soft Sweeps Review

####################
#Preparations
#setwd("~/Dropbox/SoftSweepsReview/Code/One_Locus/")
#system("./make_HIV1site_1locus")    #compile the C++ code 

listMutRates<-round(10^(seq(-7,-4,by=.5))/2,8)
#listMutRates<-listMutRates[c(2,4,6)]
options("scipen" = 999)

epistasis = 0
for (sd in c(0.001, 0.1)){
    NUMRUNS = 1000
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
        x<-c(x,paste("det=",0,sep=""))
        x<-c(x,paste("newmut=",1,sep=""))
        x<-c(x,
             "echo \"", "$seed", "$NUMRUNS", "$mu", "$sd","$sb", "$epi","$det","$newmut",
             paste("\" | ./CodeSoftSweepsReview_1locus > out_", mu, "_" , sd, "_", "1loc" ,".txt",sep="")
        )
        write(x,file="./tempscript.sh")
        system("chmod 775 ./tempscript.sh")
        #Run tempscript.sh
        #system("./tempscript.sh")
    }
}

for (sd in c(0.001, 0.1)){
    #Read the output of the sims into a dataframe called Data
    Data<-data.frame("N"=0, "Theta_trait"=0,"Theta_locus"=0,"NumHard" = 0,"NumHardSGV" = 0,"NumSos" = 0, "NumMos" = 0)
    for (i in 1:length(listMutRates)){
        #for (i in 1:2){
        mu = listMutRates[i]
        print(paste("mu",mu, "epistasis", epistasis))
        filename <- paste("out_", mu, "_" , sd, "_", "1loc" ,".txt",sep="")
        print(filename)
        x<-read.csv(filename,sep="\t",header = FALSE)
        Data[i,]<-x[1,c(2,4,6,8,10,12,14)]
    }
    
    if (sd == 0.001) filename = "Newonelocus0_001.pdf"
    if (sd == 0.1) filename = "Newonelocus0_1.pdf"
    t="A. Strong trade-off, single-locus target"
    if (sd==0.001) t= "B. Weak trade-off, single-locus target"
    #pdf(filename,width=10,height=8) #uncomment if you want to make the figure as a pdf, the legend doesn't show well in R
    par(mar = c(4,4,10,1))
    if (sd==0.1){
    barplot(rbind(Data$NumHard,Data$NumHardSGV,Data$NumSos,Data$NumMos),
            main = t, cex.main=1.8,yaxt="n",
            legend.text=c("Hard sweep - new mutation","Hard sweep - SGV", "Soft sweep - single origin", "Soft sweep - multiple origins" ),
            args.legend = list(x = 4.3, y=370,bg = "white",cex=1.5),
            col=c('#d7191c','#fdae61','#ffffbf','#abd9e9','#2c7bb6')[2:5]
    )}
    if (sd==0.001){
        barplot(rbind(Data$NumHard,Data$NumHardSGV,Data$NumSos,Data$NumMos),
                main = t, cex.main=1.8, yaxt="n",
                col=c('#d7191c','#fdae61','#ffffbf','#abd9e9','#2c7bb6')[2:5]
        )}
    mtext(c("0.001","0.003","0.01","0.03","0.1","0.3","1"), side=1,at = c(0.7,1.9,3.1,4.3,5.5,6.7,7.9),cex=1.4,line=0.5)
    mtext(expression(paste(Theta [l], " = ", Theta [g])),side=1,line=2.3,cex=1.5)
    mtext("Number of simulation runs", side=2,cex=1.5,line=2)
    axis(2,rep("",6),at=seq(0,1000,by=200),line=-1.)
    mtext(seq(0,1000,by=200),side=2,at=seq(0,1000,by=200),line=-0.5,las=1,cex=1.3)
    #dev.off() #uncomment if you want to make the figure as a pdf
}
