Getting Started
---------------

    library(VSJM)

Reading Data
-----------------------------
rm(list=ls())
Data=read.table("Data.txt",header=TRUE)
attach(Data)
names(Data)
surt=surt
event=event
Y=cbind(sbp1,sbp2,sbp3,sbp4,sbp5)
time=cbind(time1,time2,time3,time4,time5)
BMI=cbind(bmi1,bmi2,bmi3,bmi4,bmi5)
x=cbind(1,Gender,Smoking,Age,Drug2,Drug1,BMI)
w=cbind(1,Gender,Smoking,Age,Drug2,Drug1)
-------------------------

Res=DS(Y=Y, time=time, x=x, w=w, surt=surt, event=event,
       niter=1000, nburnin=500, nthin=2, nchains=2)

