library(plyr)
library(reshape2)
library(ggplot2)
library(seacarb)
library(marelac)

setwd("~/Documents/Zx‘s Papers/Arctic pCO2 science/7th round/revision based on comments/Figs/new data process 20190606/melting amplication")

#for basin
df <- data.frame(order=c(1:(107*24)),
                 year=rep(1994:2017,each=107),
                 month=(rep(rep(7:10,times=c(31,31,30,15)),times=24)),
                 day=rep(c(1:31,1:31,1:30,1:15),times=24),
                 #MLD=rep(c(15,10),times=c(1196,1012)),
                 MLD=(rep(rep(c(10,10,15,15),times=c(31,31,30,15)),times=24)),
                 NCP=(rep(rep(c(1),times=c(107)),times=24))
                 )
#for shelf
df <- data.frame(order=c(1:(107*24)),   #92 for month7-9,107 for month 7-10
                 year=rep(1994:2017,each=107),
                 month=(rep(rep(7:10,times=c(31,31,30,15)),times=24)),
                 day=rep(c(1:31,1:31,1:30,1:15),times=24),
                 day=rep(c(1:31,1:31,1:30,1:15),times=24),
                 MLD=20,
                 NCP=(rep(rep(c(5),times=c(107)),times=24))  ###for constant NCP
                 #NCP=rep(c(10,10*1.3),times=c(1391,1177))  ### for 30% increase NCP after 2006
                 )
               

simu <- read.csv("CanadaBain_seasonal simulation.csv")
simu <- read.csv("Chukchi shelf simulation.csv")
mean(simu$pCO2w)

simu$pCO2w=260 #set initial pCO2 260 for shelf 

df <- merge(df,simu,by=c("year","month"),all.x=TRUE)
df <- df[order(df$order),]
df$date <- df$year+(df$month-0.5)/12

####interpolation of ice melting
ice <- data.frame()
icey <- unique(df$year)
for(i in 1:length(icey)){
        I <- data.frame(ice=c(1:107))
        icea <- subset(simu,simu$year==icey[i])
        I$ice[1:15]<- icea$ICE...[1]
        t1 <- c(1:31)
        dat1 <- data.frame(x=c(1,31),y=icea$ICE...[1:2])
        dat2 <- approx   (dat1$x, dat1$y , xout=t1, method = "linear", n = 31)
        I$ice[16:46] <- dat2$y
        t2 <- c(1:31)
        dat1 <- data.frame(x=c(1,31),y=icea$ICE...[2:3])
        dat2 <- approx   (dat1$x, dat1$y , xout=t2, method = "linear", n = 31)
        I$ice[47:77] <- dat2$y
        t3 <- c(1:30)
        dat1 <- data.frame(x=c(1,30),y=icea$ICE...[3:4])
        dat2 <- approx   (dat1$x, dat1$y , xout=t3, method = "linear", n = 30)
        I$ice[78:107] <- dat2$y
        ice <- rbind(ice,I)
}

plot(ice$ice)
df$ICE... <- ice$ice

#df$ICE... <- 0

###SST smoothing
temp <- data.frame()
icey <- unique(df$year)
for(i in 1:24){
        T <- data.frame(temp=c(1:107))
        Ta <- subset(simu,simu$year==icey[i])
        T$temp[1:15]<- Ta$SST.degC. [1]
        t1 <- c(1:31)
        dat1 <- data.frame(x=c(1,31),y=Ta$SST.degC.[1:2])
        dat2 <- approx   (dat1$x, dat1$y , xout=t1, method = "linear", n = 31)
        T$temp[16:46] <- dat2$y
        t2 <- c(1:31)
        dat1 <- data.frame(x=c(1,31),y=Ta$SST.degC.[2:3])
        dat2 <- approx   (dat1$x, dat1$y , xout=t2, method = "linear", n = 31)
        T$temp[47:77] <- dat2$y
        t3 <- c(1:30)
        dat1 <- data.frame(x=c(1,31),y=Ta$SST.degC.[3:4])
        dat2 <- approx   (dat1$x, dat1$y , xout=t3, method = "linear", n = 30)
        T$temp[78:107] <- dat2$y
        temp <- rbind(temp,T)
}

mean(df$SST.degC.)

df$SST.degC. <- temp$temp+0.30 ##+ 0.63 for basin and + 0.30 for shelf

#df$SST.degC.=mean(df$SST.degC.)+0.63  #-0.18 for Basin
plot(df$SST.degC.)
summary(lm(df$SST.degC.~df$date))
#sal derived from Ice
#for(i in 1:length(df$ICE...)){
#        if( df$ICE...[i]<80) {
#                df$Sal[i] <- 3.75*df$ICE...[i]/100+26
#        } else {
#                df$Sal[i] <- 5*df$ICE...[i]/100+25
#        }     
#        }

#df$Sal=27.3 #for bain
#df$Sal=0.039*df$ICE...+25.6  #for basin
#df$TA=df$Sal*56.20+432# for Basin from GLODAP 2019 data
#df$DIC=1865/27.2*df$Sal

df$Sal=29.8 #for shelf
#df$Sal=0.027*df$ICE...+29.5 #for shelf
df$TA=df$Sal*59.41+310 ###for shelf 

plot(df$date, df$Sal)
mean(df$Sal)
sd(df$Sal)
summary(lm(df$Sal~df$date))
#simu_monthlymean <- aggregate(df,mean,na.rm=TRUE, by=list(month=df$month,year=df$year))
#simu_monthlymean<- simu_monthlymean[,-c(1:2)]
#simu_monthlymean$DATE <- simu_monthlymean$year+(simu_monthlymean$month-0.5)/12

#p <- ggplot()
#p+geom_point(data=simu_monthlymean,aes(simu_monthlymean$DATE,simu_monthlymean$Sal))+
#        geom_point(data=monthlymean,aes(monthlymean$Date_year,monthlymean$SAL),color="red")






a<- carb(flag=24,var1=df$pCO2w,var2=df$TA/1000000,S=df$Sal,T=df$SST.degC., Patm=1, P=0, Pt=0, Sit=0,
         k1k2="l", kf="dg", ks="d", pHscale="T", b="l10", gas="insitu", warn="y")
df$DIC=a$DIC*1000000

df$sDIC=df$DIC/df$Sal*29.8 #### 27.3 for basin ; 29.8 for shelf
plot(df$sDIC)
summary(lm(df$sDIC~df$date))
#a<- carb(flag=15,var1=df$TA/1000000,var2=df$DIC/1000000,S=df$Sal,T=df$SST.degC., Patm=1, P=0, Pt=0, Sit=0,
#         k1k2="l", kf="dg", ks="d", pHscale="T", b="l10", gas="insitu", warn="y")
#df$pCO2w=a$pCO2insitu


plot(df$Sal)
plot(df$TA)
plot(df$DIC)
plot(df$pCO2w)
#######calculation of Sc number
A=2116.8
B=-136.25
C=4.7353
D=-0.092307
E=0.0007555
##Sc=A+B*t+C*t^2+D*t^3+E*t^4 (t in C)
##t=20 Sc=668
df$Sc <- A+B*df$SST.degC.+C*df$SST.degC.^2+D*df$SST.degC.^3+E*df$SST.degC.^4

###calculate the solubility of CO2 unit in moles/l atm
A1=-58.0931
A2=90.5069
A3=22.2940
B1=0.027766
B2=-0.025888
B3=0.0050578
## ln(b)=A1+A2*(100/T)+A3*log(T/100)+S*(B1+B2*(T/100)+B3*(T/100)^2)
Bunsen=A1+A2*(100/(df$SST.degC.+273.15))+A3*log((df$SST.degC.+273.15)/100)+df$Sal*(B1+B2*((df$SST.degC.+273.15)/100)+B3*((df$SST.degC.+273.15)/100)^2)
df$solubility <- exp(Bunsen)
df$density <- rho(S=df$Sal,T=df$SST.degC.,P=0)
df$solubility <- df$solubility/(df$density/1000)


df2 <- subset(df,df$month==9)
plot(df2$SST.degC.)


####calculate the gas tranfer volocity of CO2
df$K <- 0.251*df$WindSpeed.m.s.*(df$Sc/660)^(-0.5)*((100-df$ICE...)/100)

df$flux <- df$K*df$solubility*(df$pCO2w-df$pCO2air)*24/100

df$deltaDIC=df$flux/df$MLD##### includes both air-sea exchange and NCP alteration
df$deltaDICb=df$NCP/df$MLD 
plot(df$deltaDIC)



Year=c(1994:2017)
df1<- data.frame() 
df2 <- data.frame() 
for( j in 1:length(Year)){
        df1 <- subset(df,df$year==Year[j])
        df1$DRDIC[1]=0
        df1$dTA[1]=0
for( i in 1:(length(df1$year)-1)) {
        
        DR <- (1+(df1$TA[i+1]-df1$TA[i])/df1$TA[i])
        df1$DIC[i+1]=df1$DIC[i]*DR  -df1$deltaDIC[i]-df1$deltaDICb[i]
        df1$DRDIC[i+1]=df1$DIC[i]*DR-df1$DIC[i]
        df1$dTA[i+1]=df1$TA[i+1]-df1$TA[i]
        b<- carb(flag=15,var1=(df1$TA[i+1]/1000000),var2=(df1$DIC[i+1]/1000000),S=df1$Sal[i+1],T=df1$SST.degC.[i+1], Patm=1, P=0, Pt=0, Sit=0,
                 k1k2="l", kf="dg", ks="d", pHscale="T", b="l10", gas="insitu", warn="y")
        
        df1$pCO2w[i+1] <- b$pCO2
        df1$flux[i+1]=df1$K[i+1]*df1$solubility[i+1]*(df1$pCO2w[i+1]-df1$pCO2air[i+1])*24/100
        df1$deltaDIC[i+1]=df1$flux[i+1]/df1$MLD[i+1]
        #df1$pH[i+1]=b$pH
        #df1$Omega[i+1]=b$OmegaAragonite
}

     
df2 <- rbind(df2,df1) 

}

#####
plot(df2$pCO2w)
plot(df2$Sal)

#plot(df2$DIC)
#plot(df2$DIC/df2$Sal*27.3-mean(df2$DIC/df2$Sal*27.3)) # for basin
#df2$DsDIC <- df2$DIC/df2$Sal*27.3-mean(df2$DIC/df2$Sal*27.3)

plot(df2$DIC/df2$Sal*29.8)  # for shelf
plot(df2$DIC/df2$Sal*29.8-mean(df2$DIC/df2$Sal*29.8)) # for shelf
df2$DsDIC <- df2$DIC/df2$Sal*29.8-mean(df2$DIC/df2$Sal*29.8)

df2$DATE=paste(df2$month,df2$day,df2$year,sep="/")
df2$DATE <- as.Date(df2$DATE,format="%m/%d/%Y")
df2$dd <- df2$year+(df2$month-0.5)/12

df3 <- aggregate(df2,mean,by=list(df2$month,df2$year))
df3$sDIC <- df3$DIC/df3$Sal*mean(df3$Sal)
plot(df3$year,df3$deltaDIC)
plot(df3$year,df3$pCO2w)

df3$day=15
df3$DATE=paste(df3$month,df3$day,df3$year,sep="/")
df3$DATE <-  as.Date(df3$DATE,format="%m/%d/%Y")
model <- lm(df3$pCO2w~df3$DATE)
model$coefficients[2]*365.25


df3$dd <- df3$year+(df3$month-0.5)/12
summary(lm(df3$pCO2w~df3$dd))
summary(lm(df3$sDIC~df3$dd))
summary(lm(df3$DsDIC~df3$dd))
df3 <- df3[,-c(1:2)]
#df3 <- subset(df3,df3$year%in%c(1994,1997,1999:2000,2003:2004,2008:2012,2015:2017)) #for Basin

#df3 <- subset(df3,df3$year%in%c(1994,1998:2000,2003:2004,2006,2008:2017)) #for shelf

#df3 <- subset(df3,!df3$month==c(7))
summary(lm(df3$pCO2w~df3$dd))


####plot + observaed monthlymean
setwd("~/Documents/Zx‘s Papers/Arctic pCO2 science/7th round/revision based on comments/Figs/new data process 20190606")
dfall<- read.csv("gridbasin(0.25*0.1)_20190606.csv")
dfall<- read.csv("gridshelf(0.25*0.1)_20190606.csv")
monthlymean <- aggregate(dfall,mean,na.rm=TRUE, by=list(month=dfall$month,year=dfall$year))
summary(lm(monthlymean$pCO2_sst~monthlymean$Date_year))

p <- ggplot()

p+geom_point(data=df2,aes(x=df2$dd,df2$pCO2w),color="grey",size=0.3)+
        geom_point(data=df3,aes(x=df3$dd,df3$pCO2w),color="black", size=2)+
        geom_smooth(data=df3,aes(x=df3$dd,df3$pCO2w),method="lm",linetype="solid",color="black",se=F,size=0.7)+
        geom_point(data=monthlymean,aes(x=monthlymean$Date_year,monthlymean$pCO2_sst),color="red",size=1)+
        #geom_point(data=df3n,aes(x=df3n$DATE,df3n$pCO2w),color="red", size=2)+
        #geom_smooth(data=df3n,aes(x=df3n$DATE,df3n$pCO2w),method="lm",linetype="solid",color="red",se=F)+
        #scale_y_continuous(breaks = seq(200,400, by=50),limits=c(200,400)) +   # for shelf 190,350  for basin 
        scale_x_continuous(breaks = seq(1994,2018, by=1),limits=c(1994,2018)) +   # for basin 200,400
        scale_y_continuous(breaks = seq(100,350, by=50),limits=c(120,320)) +   # for shelf 190,350  for shelf
        theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

p <- ggplot()

p+#geom_point(data=df2,aes(x=df2$DATE,df2$sDIC),color="grey",size=0.5)+
        geom_point(data=df3,aes(x=df3$dd,df3$DsDIC),color="black", size=2)+
        geom_smooth(data=df3,aes(x=df3$dd,df3$DsDIC),method="lm",linetype="solid",color="black",se=F,size=0.5)+
        #geom_point(data=df2n,aes(x=df2n$DATE,df2n$pCO2w),color="grey",size=0.3)+
        #geom_point(data=df3n,aes(x=df3n$DATE,df3n$pCO2w),color="red", size=2)+
        #geom_smooth(data=df3n,aes(x=df3n$DATE,df3n$pCO2w),method="lm",linetype="solid",color="red",se=F)+scale_y_continuous(breaks = seq(200,400, by=50),limits=c(200,400)) +   # for shelf 190,350  for basin 
        scale_x_continuous(breaks = seq(1994,2018, by=1),limits=c(1994,2018)) +   # for basin 200,400
        
        theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))




        



##500 X 400



###last one -first one
d <- data.frame()
for( k in 1:length(Year)){
        df4 <- subset(df2,df2$year==Year[k])
        dpCO2 <- df4$pCO2w[92]-df4$pCO2w[1]
        
        d <- rbind(d,dpCO2)
        
}

colnames(d)<-c("dpCO2")

##max -min choose this one
d <- data.frame()
for( k in 1:length(Year)){
        df4 <- subset(df2,df2$year==Year[k])
        dpCO2 <- df4$pCO2w[which.max(df4$pCO2w)] - df4$pCO2w[which.min(df4$pCO2w)]
        
        d <- rbind(d,dpCO2)
        
}

colnames(d)<-c("dpCO2")



df3 <- subset(df3,df3$month==9)


df4 <- data.frame(year=c(1994:2017),dpCO2=d$dpCO2)
df4$DATE=df3$DATE
df4$ice=df3$ICE...
plot(df4$year,df4$dpCO2)
lm(df4$dpCO2~df4$year)
anova(lm(df4$dpCO2~df4$year))
summary(lm(df4$dpCO2~df4$year))
summary(anova(lm(df4$dpCO2~df4$year)))

lm(df4$ice~df4$dpCO2)
summary(lm(df4$ice~df4$dpCO2))
plot(df4$ice,df4$dpCO2)
p <- ggplot(data=df4)
p+geom_point(aes(x=df4$ice,y=df4$dpCO2,color=df4$year),size=4)+
        scale_colour_gradient(name="Year",low = "grey", high = "orange")+
        #geom_smooth(aes(x=df4$ice,y=df4$dpCO2),method="lm")+
        #geom_line(aes(x=df4$ice,y=df4$dpCO2),size=0.7)+
        xlab("Ice%")+
        ylab("Seasonal increased pCO2")+
        #scale_x_continuous(breaks = seq(1994,2018, by=2),limits=c(1994,2018)) + 
        #scale_y_continuous(breaks = c(150,200,250,300,350,400),limits=c(260,405))+
        #geom_vline(xintercept=c(10,79,130), linetype="dashed",color="darkblue")+
        #geom_hline(yintercept=c(378,404), linetype="dashed",color=c("blue","darkblue"),size=1)+
        #scale_x_continuous(breaks = c(0,30,61,93,123,153),limits=c(0,153))+
        theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#600 X 400
#800x500


#df4$year <- as.numeric(df4$year)
#df4$dpCO2 <- as.numeric(df4$dpCO2)


p <- ggplot()
p+geom_point(data=df4,aes(x=df4$year,y=df4$dpCO2,color=df4$ice),size=4)+
        geom_smooth(data=df4,aes(x=df4$year,y=df4$dpCO2),method='lm',color="orange",se=F)+
        scale_colour_gradient(name="Ice%",low = "grey", high = "darkblue")+
        #geom_errorbar(aes(x=df3$DATE,ymin=(df3$ICE...-df3$ICE_stdRg3), ymax=(df3$ICE...+df3$ICE_stdRg3),width=3))+
        geom_line(aes(x=df4$year,y=df4$dpCO2),size=0.7,color="orange")+
        xlab("year")+
        ylab("Seasonal increased pCO2")+
       scale_x_continuous(breaks = seq(1994,2018, by=2),limits=c(1994,2018)) + 
        scale_y_continuous(breaks = c(0,25,50,75),limits=c(0,76.5))+
        #geom_vline(xintercept=c(10,79,130), linetype="dashed",color="darkblue")+
        #geom_hline(yintercept=c(378,404), linetype="dashed",color=c("blue","darkblue"),size=1)+
        #scale_x_continuous(breaks = c(0,30,61,93,123,153),limits=c(0,153))+
        theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
        geom_point(data=dfNCP10,aes(x=dfNCP10$year,y=dfNCP10$dpCO2,color=dfNCP10$ice),size=4)+
        geom_smooth(data=dfNCP10,aes(x=dfNCP10$year,y=dfNCP10$dpCO2),method='lm',color="black",se=F)+
        geom_line(aes(x=dfNCP10$year,y=dfNCP10$dpCO2),size=0.7,color="black")
        
        

#Size= 800 x 500

dfNCP10 <- df4
dfNCP12 <- df4


###for plotting
df3$DATE=paste(df3$month,df3$day,df3$year,sep="/")
df3$DATE <-  as.Date(df3$DATE,format="%m/%d/%Y")

model <- lm(df3$pCO2w~df3$DATE)
model$coefficients[2]*365.25




plot(df3$DATE,df3$pCO2w)


d <- data.frame()
for( k in 1:length(Year)){
        df5 <- subset(df2,df2$year==Year[k])
        dDIC <- df5$DIC[107]-df5$DIC[1]
        
        d <- rbind(d,dDIC)
        
}

colnames(d)<-c("dDIC")

df2 <- df2[,-c(33)]
df5 <- aggregate(df2,sum,by=list(df2$year))
df5$DATE=df4$DATE
df6 <- data.frame(
                  year=df5$Group.1,
                  DICairflux=-df5$deltaDIC,
                  DRDIC=df5$DRDIC,
                  DICb=-df5$deltaDICb
                )

df7 <- melt(df6,id = c("year"))


df6$netDIC=df6$DRDIC+df6$DICairflux+df6$DICb
df6$netTA=df5$dTA
df6$dDICdTA=df6$netDIC/df6$netTA
plot(df6$year,df6$dDICdTA)
#df6$year <- as.character(df6$year)


p <- ggplot()
p+geom_bar(data=df7,aes(x=df7$year,y=df7$value,fill=df7$variable),stat="identity")+ #
        geom_line(data=df6,aes(x=df6$year,y=df6$netDIC),color="brown2")+
        geom_point(data=df6,aes(x=df6$year,y=df6$netDIC),size=2,color="brown2")+
        geom_point(data=df6,aes(x=df6$year,y=df6$netTA),size=2,color="black")+
        geom_line(data=df6,aes(x=df6$year,y=df6$netTA),color="black")+
        #geom_errorbar(aes(x=df3$DATE,ymin=(df3$ICE...-df3$ICE_stdRg3), ymax=(df3$ICE...+df3$ICE_stdRg3),width=3))+
        #geom_line(aes(x=df5$DATE,y=df5$deltaDIC),size=0.7)+xlab("day")+
        #geom_smooth(aes(x=df5$DATE,y=df5$deltaDIC),method=lm,se=FALSE,size=1,color="darkblue")+ ##normal size=2
        ylab("uptake DIC)")+
        scale_x_continuous(breaks = seq(1994,2018, by=2),limits=c(1993.5,2018)) + 
        #scale_y_continuous(breaks = c(150,200,250,300,350,400),limits=c(260,405))+
        #geom_vline(xintercept=c(10,79,130), linetype="dashed",color="darkblue")+
        #geom_hline(yintercept=c(378,404), linetype="dashed",color=c("blue","darkblue"),size=1)+
        #scale_x_continuous(breaks = c(0,30,61,93,123,153),limits=c(0,153))+
        scale_fill_brewer(palette="Paired")+
        theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))












#size 800 X 500
#for Chukchi Sea 480625 km2
#for Canada basin 983750 km2

df5$uptake <- 480625*1000000*df5$flux/1000*12/(10^12)
df5$uptake <- 983750*1000000*df5$flux/1000*12/(10^12)

plot(df5$uptake)


df8 <- data.frame(year=df5$Group.1,uptakeh=-df5$uptake)
df9 <- data.frame(year=c(1998:2000,2003,2004,2006,2008:2017),
                  uptakel=c(3.77,4.68,3.83,6.81,5.23,5.73,6.19,5.17,8.25,3.20,7.98,5.91,8.66,6.46,7.68,7.39))
df9 <- data.frame(year=c(1998:2000,2003,2004,2006,2008:2017),
                  uptakel=c(2.38,2.38,1.33,3.91,3.05,3.05,5.08, 3.54,4.31, 3.68, 4.50, 1.48, 2.83,3.98,4.14, 4.85))



df8 <- merge(df8,df9,by=c("year"),all.x=TRUE)

p <- ggplot(data=df8)
p+geom_ribbon(aes(ymin=df8$uptakel,ymax=df8$uptakeh,x=df8$year,fill="band"),alpha=0.5)+
        geom_line(aes(x=df8$year,y=df8$uptakeh))+
        geom_line(aes(x=df8$year,y=df8$uptakel))
        


