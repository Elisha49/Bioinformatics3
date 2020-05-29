#Qsn no:1
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part1_files/gene_expression.tsv",
     destfile = "gene_expression.tsv")

x<-read.table("gene_expression.tsv")

head(x)
str(x)
x <- read.csv("gene_expression.tsv", sep="\t",stringsAsFactors = FALSE)
head(x)
str(x)
x <- read.table("gene_expression.tsv",header = TRUE , stringsAsFactors = FALSE , row.names = 1)
head(x)
str(x)
#Qsn no: 2
x$meana <- rowMeans(x)
head(x)
#Qsn no :3
order(x$meana)
order(-x$meana)
x[order(-x$meana),]
x_sorted<-x[order(-x$meana),]
head(x_sorted,10)
#Qsn no : 4
head(x,10)
which(x[,3]>10)
dim((x[which(x[,3]<10),]))
#QNA 5
hist(x$meana)
hist(x$meana,breaks = 10)

