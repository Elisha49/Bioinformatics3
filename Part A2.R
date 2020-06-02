#Qsn no:6(To download the file "growth_data.csv")
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part1_files/growth_data.csv",
              destfile ="growth_data.csv")
# To read the file in R
abc<-read.table("growth_data.csv",sep = ",",header = TRUE,stringsAsFactors = FALSE)
head(abc)
str(abc)
# to know the colnames of the file
colnames(abc)

#Qsn no: 7(subset function of the sites)
northeast <- subset(abc, Site =="northeast")
head(northeast)
tail(northeast)
southwest <- subset(abc, Site =="southwest")
head(southwest)
tail(southwest)
# To calculate the mean and standard deviation of the sites
mean(northeast$Circumf_2004_cm)
mean(northeast$Circumf_2019_cm)
sd(northeast$Circumf_2004_cm)
sd(northeast$Circumf_2019_cm)
mean(southwest$Circumf_2004_cm)
mean(southwest$Circumf_2019_cm)
sd(southwest$Circumf_2004_cm)
sd(southwest$Circumf_2019_cm)

#Qsn no:8 
northeast <- subset(abc, Site=="northeast")
southwest <- subset(abc, Site=="southwest")
# To make a box plot
boxplot(northeast$Circumf_2004_cm)
boxplot(northeast$Circumf_2004_cm[1:50],ylab="Circumference(cm)",main="Growth at northeast 2019")
boxplot(northeast$Circumf_2019_cm)
boxplot(northeast$Circumf_2019_cm[1:50],ylab="Circumference(cm)",main="Growth at northeast 2004")
boxplot(southwest$Circumf_2004_cm)
boxplot(southwest$Circumf_2004_cm[1:50],ylab="Circumference(cm)",main="Growth at southwest 2004")
boxplot(southwest$Circumf_2019_cm)
boxplot(southwest$Circumf_2019_cm[1:50],ylab="Circumference(cm)",main="Growth at southwest 2019")
#Boxplot of sites
boxplot(northeast$Circumf_2004_cm,northeast$Circumf_2019_cm,southwest$Circumf_2004_cm,southwest$Circumf_2019_cm)
#Boxplot of sites with cm and main
boxplot(northeast$Circumf_2004_cm,northeast$Circumf_2019_cm,southwest$Circumf_2004_cm,southwest$Circumf_2019_cm,
        names = c("northeast2004","northeast2019","southwest2004","southwest2019"),
        ylab="Circumference(cm)",main="Growth at start and end year")
#Qsn no:9
#(To calculate the mean growth of northeast site at Circumf_2009_cm) 
mean(northeast$Circumf_2009_cm)
x <- mean(northeast$Circumf_2009_cm)
#(To calculate the mean growth of northeast site at Circumf_2019_cm)
mean(northeast$Circumf_2019_cm)
y <- mean(northeast$Circumf_2019_cm)
a <- sum(x,y) 
mean <- a/2
str(northeast)

#(To calculate the mean growth of southwest site of circumf_2009_cm )
mean(southwest$Circumf_2009_cm)
n <- (southwest$Circumf_2009_cm)
#(To calculate the mean growth of southwest site of circumf_2019_cm )
mean(southwest$Circumf_2019_cm)
m <- mean(southwest$Circumf_2019_cm)
b <- sum(n,m)
mean <- b/2
#(To calculate the mean growth of northeast site over the past 10 years )
northeastmean <- northeast$Circumf_2019_cm-northeast$Circumf_2009_cm
northeast$growth <- northeast$Circumf_2019_cm-northeast$Circumf_2009_cm
mean(northeast$growth)
head(northeast$growth)
head(northeast)
str(northeast)
#(To calculate the mean growth of southwest site over the past 10 years )
southwestmean <- southwest$Circumf_2019_cm-southwest$Circumf_2009_cm
southwest$growth <- southwest$Circumf_2019_cm-southwest$Circumf_2009_cm
mean(southwest$growth)
head(southwest$growth)
head(southwest)
str(southwest)
#Qsn no:10(Using t.test and wilcox.test to estimate the p-value)
t.test(northeast$growth,southwest$growth)
wilcox.test(northeast$growth,southwest$growth)
