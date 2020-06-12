#Qsn no:1
#To download the file "gene_expression.tsv"
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part1_files/gene_expression.tsv",
              destfile = "gene_expression.tsv")

# to read the file in R
x<-read.table("gene_expression.tsv")
head(x)
str(x)
x <- read.csv("gene_expression.tsv", sep="\t",stringsAsFactors = FALSE)
head(x)
str(x)

# to show the table for first six genes
x <- read.table("gene_expression.tsv",header = TRUE , stringsAsFactors = FALSE , row.names = 1)
head(x)
str(x)

#Qsn no:2(to make a new coloumn which is the mean of other column)
x$meana <- rowMeans(x)
head(x)

#Qsn no:3(To list the genes with highest mean expression)
order(x$meana)
order(-x$meana)
x[order(-x$meana),]
#To list the top 10 genes with highest mean expression
x<-x[order(-x$meana),]
head(x,10)

#Qsn no:4(To determines the number of genes with mean<10)
y <- subset(x, meana<10) 
# Number of genes with mean<10
head(y)
str(y)

#Qsn no:5(histogram plot of the mean values)
hist(x$meana)
hist(x$meana,breaks = 20)

