########### TABLE OF CONTENTS ############################
# 0.  Define all libraries and functions
# 1.  Human inputs - manually enter file name as variable f (USER INPUT: file name)
# 2.  File 
 
#0. Define all libraries used and functions
# None required yet


# Define functions for converting oxygen
o2sol <- function(T=0, S=0, Pb=1013)
{
  # derived from eqn 8 in Garcia and Gordon, 1992, Limnol Oceanog
  # default values are 0oC and 0 salinity
  # solubility values are given for air-saturated water, at STP (760 mmHg
  # total atmospheric temperature, including WVP)
  # solubility units are in umol O2/kg water
  # T=water temperature (°C), S=salinity (%o, per mil, ppt)
  # Pb=barometric pressure in mbar.  Solubility will be corrected to this pressure condition
  A0 <- 5.80818
  A1 <- 3.20684
  A2 <- 4.11890
  A3 <- 4.93845
  A4 <- 1.01567
  A5 <- 1.41575
  B0 <- -7.01211e-3
  B1 <- -7.25958e-3
  B2 <- -7.93334e-3
  B3 <- -5.54491e-3
  C0 <- -1.32412e-7
  Ts<-log((298.15-T)/(273.15+T))
  ln.o2sol<-A0 + A1*Ts + A2*Ts^2 + A3*Ts^2 + A3*Ts^3 + A4*Ts^4 + A5*Ts^5 + S*(B0 + B1*Ts + B2*Ts^2 + B3*Ts^3) + C0*S^2
  o2sol<-exp(ln.o2sol)*Pb/1013
  o2sol
}

pH2O<-function(T=20)
{
  # default temperature is 20°C
  # pH2O units in mbar
  pH2O <- 6.112 * exp ( 17.62*T / (243.12 + T))
  pH2O
}

#1. Human Inputs - Filename, respirometer volumes, animal mass ####
mainDir<-"~/Documents/R/MyProjects/PyroscienceImport/data"
setwd(mainDir)
getwd() #show working directory
l.files<-list.files()
# define the file output names and locations
#  must give a filename ("in_file")

f<-"low flow test may 29.txt"
f<-"la cockaroochasmooth.txt"

vol.units<-"mL"
mass.units<-"mg"
# input respirometer volumes here, including tubing and pumps
v1<-1
v2<-1
v3<-1
v4<-1
# input animal masses here
m1<-1
m2<-1
m3<-1
m4<-1
  

# Header Info
head<-read.csv(f, header=FALSE, sep="", nrow=20) 
smoothed<-as.numeric(as.character(head[16,3]))
smoothed # no. of samples data were smoothed


# Read alldata information
alldata<-read.delim(f, header=TRUE, sep="\t", skip=19) 
head(alldata)
colnames(alldata)<-NULL
colnames(alldata)<-c(alldata[1,])
date.time<-strptime(paste(alldata[,1],alldata[,2]), format="%d.%m.%Y %H:%M:%S", tz="EST5EDT")
alldata<-alldata[,-c(1,2)]
alldata<-data.frame(date.time,alldata)
alldata<-alldata[-c(14:ncol(alldata))]
# remove all the extra columns pyroscience generates.  

colnames(alldata)<-c("Date.Time","Elapsed Time (s)", "Comments", "Ch1-O2 (µmol/L)", "Ch2-O2  (µmol/L)","Ch3-O2  (µmol/L)",
                     "Ch4-O2  (µmol/L)", "Ch1-T  (°C)", "Ch2-T  (°C)", "Ch3-T  (°C)", "Ch4-T  (°C)", "Pressure (mbar)", 
                     "Humidity (%)") 
#alldata<-na.omit(alldata)

# Smooth the data using the lowess function
smoother.span<-1/20
smooth.ch1<-lowess(alldata[,2],alldata[,4],f=smoother.span)$y
smooth.ch2<-lowess(alldata[,2],alldata[,5],f=smoother.span)$y
smooth.ch3<-lowess(alldata[,2],alldata[,6],f=smoother.span)$y
smooth.ch4<-lowess(alldata[,2],alldata[,7],f=smoother.span)$y

# plot the raw and smooth data together
par(mfrow=(c(2,2)), mar=c(4,4,2,2), cex=0.8)
plot(alldata[,2],alldata[,4], 
     xlab=colnames(alldata)[2],
     ylab=colnames(alldata)[4])
points(alldata[,2],smooth.ch1,type="l",col="red",lwd=4)
plot(alldata[,2],alldata[,5], 
     xlab=colnames(alldata)[2],
     ylab=colnames(alldata)[5])
points(alldata[,2],smooth.ch2,type="l",col="red",lwd=4)
plot(alldata[,2],alldata[,6], 
     xlab=colnames(alldata)[2],
     ylab=colnames(alldata)[6])
points(alldata[,2],smooth.ch3,type="l",col="red",lwd=4)
plot(alldata[,2],alldata[,7], 
     xlab=colnames(alldata)[2],
     ylab=colnames(alldata)[7])
points(alldata[,2],smooth.ch4,type="l",col="red",lwd=4)

# calculate slopes (mult by -1) of the smoothed data
slope.ch1<- -1*diff(smooth.ch1)/diff(alldata[,2])
slope.ch2<- -1*diff(smooth.ch2)/diff(alldata[,2])
slope.ch3<- -1*diff(smooth.ch3)/diff(alldata[,2])
slope.ch4<- -1*diff(smooth.ch4)/diff(alldata[,2])
diff.time<-alldata[-1,2]


# Significant Events - Denoted by non-zero comments
comments.time<-alldata[which(!alldata[,3]==""),2]
comments<-as.character(alldata[which(!alldata[,3]==""),3])


# Plot the final data 
par(mfrow=(c(2,2)), mar=c(5,5,3,2), cex=0.8)
plot(diff.time, slope.ch1, type="l", main="Channel 1 (µmol/L/s)",
     xlab=colnames(alldata)[2], ylab="Slope (µmol/L/s)", bty="n")
rect(0,-1,max(diff.time),0,col="light grey", angle=90, border=1)
lines(diff.time, slope.ch1)
text(comments.time, max(slope.ch1), substr(comments,1,1), col="red")

plot(diff.time, slope.ch2, type="l", main="Channel 2 (µmol/L/s)",
     xlab=colnames(alldata)[2], ylab="Slope", bty="n")
rect(0,-1,max(diff.time),0,col="light grey", angle=90, border=1)
lines(diff.time, slope.ch2)
text(comments.time, max(slope.ch2), substr(comments,1,1), col="red")

plot(diff.time, slope.ch3, type="l", main="Channel 3 (µmol/L/s)",
     xlab=colnames(alldata)[2], ylab="Slope", bty="n")
rect(0,-1,max(diff.time),0,col="light grey", angle=90, border=1)
lines(diff.time, slope.ch3)
text(comments.time, max(slope.ch3), substr(comments,1,1), col="red")

plot(diff.time, slope.ch4, type="l", main="Channel 4 (µmol/L/s)",
     xlab=colnames(alldata)[2], ylab="Slope", bty="n")
rect(0,-1,max(diff.time),0,col="light grey", angle=90, border=1)
lines(diff.time, slope.ch4)
text(comments.time, max(slope.ch4), substr(comments,1,1), col="red")

