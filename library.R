library("seqinr")
library("R.utils")
library("rBLAST")
library("ape")
library("ORFik")
library("Biostrings")
#To download the e.coli file
download.file("ftp://ftp.ensemblgenomes.org/pub/bacteria/release-42/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/cds/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz",
              destfile = "e.coli.fa.gz")
R.utils::gunzip("e.coli.fa.gz", overwrite= TRUE)
makeblastdb("e.coli.fa",dbtype="nucl","-parse_seqids")
#To download sample file
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/sample.fa",
              destfile = "sample.fa")
x <- read.fasta("sample.fa")
myseq <- x[[2]]
myseq
str(myseq)
seqinr::getLength(myseq)
length(myseq)
seqinr::GC(myseq)
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/mutblast_functions.R",
                  destfile = "mutblast.R")

source("mutblast.R")



myblastn_tab <- function(myseq,db) {
  mytmpfile1<-tempfile()
  mytmpfile2<-tempfile()
  write.fasta(myseq,names=attr(x,"name"),file.out = mytmpfile1)
  system2(command = "/usr/bin/blastn",
          args = paste("-db ", db ," -query", mytmpfile1 ,"-outfmt 6 -evalue 0.05 -ungapped >" 
                       , mytmpfile2))
  # probably need to add headers here
  res <- NULL
  if (file.info(mytmpfile2)$size > 0 ) {
    res <- read.csv(mytmpfile2,sep="\t",header=FALSE)
    colnames(res) <- c("qseqid","sseqid","pident","length","mismatch","gapopen",
                       "qstart","qend","sstart","send","evalue","bitscore")
  }
  unlink(c(mytmpfile1,mytmpfile2))
  if (!is.null(res)  ) {
    res <- res[order(-res$bitscore),]
  }
  res
}

source("mutblast.R")
res <- myblastn_tab(myseq = x, db="e.coli.fa")
str(res)
res
head(res)

hits <- as.character(res$sseqid[1:3])
hits
db <- read.fasta("e.coli.fa")
str(db[1:6])
head(names(db))

myseqs <- db[which(names(db)%in%hits)]
myseqs <- c(myseqs,x)
seqinr::write.fasta(tophit,names = names(myseq),file.out = "myseqs.fa")
str(myseqs)

tophit <- db[which(names(db) %in% hits[1])]
tophit[1:3]
seqinr::write.fasta(tophit,names=names(tophit),file.out = "tophit.fa")
str(tophit)

makeblastdb("tophit.fa",dbtype="nucl", "-parse_seqids" )
res <- myblastn(myseq = x, db = "tophit.fa" )
res


mutator <- function(myseq,nmut) {
  myseq_mod <- myseq
  mypos<-sample(seq_along(myseq),nmut)
  myseq_mod[mypos] <- sample(c("a","c","g","t"),length(mypos),replace = TRUE)
  return(myseq_mod)
}

#qsn 4
#qsn 4
source("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/mutblast_functions.R")
# mutator function  to replace DNA bases with a random base

mutator <- function(myseq,nmut) {
  myseq_mod <- myseq
  mypos<-sample(seq_along(myseq),nmut)
  myseq_mod[mypos] <- sample(c("a","c","g","t"),length(mypos),replace = TRUE)
  return(myseq_mod)
}

x_mut <-mutator(myseq = myseq, 100)

x_mut

myblastn_tab(myseq=myseq,db="tophit.fa")
res <- myblastn_tab(myseq=myseq,db="tophit.fa")
head(res)

#Pairwise alignment created
x_mut_1 <- DNAString(c2s(x_mut))

myseq_1 <- DNAString(c2s(myseq))

aln <- pairwiseAlignment(myseq_1,x_mut_1)


pid(aln)
# To know the number of mismatches
nmismatch(aln)

#Qsn no:- 5
# Function to randomise the sequence
myfunc <- function(myseq,nmut) {
  mutseq <- mutator( myseq= myseq, nmut = nmut) 
  res <- myblastn_tab(myseq= mutseq, db= "tophit.fa") #for blast
  if (is.null(res)) {myres= 0} else {myres = 1}
  return(myres)
}

myfunc(myseq,525)
# TO replicate the sequence

replicate(50, myfunc(myseq,100) )

replicate(70, myfunc(myseq,100) )

replicate(80, myfunc(myseq,100) )

replicate(90, myfunc(myseq,100) )

replicate(100, myfunc(myseq,100) )

replicate(20, myfunc(myseq,200) )

replicate(50, myfunc(myseq,200) )

mean(replicate(50, myfunc(myseq,200) ) )

#Qsn no:- 6
nmut <- 200
variables <- 2
output <- matrix(nrow = nmut, ncol = variables)

{
  propmatch <- mean(replicate(50, routine(myseq, 200)))
  output[i, ] <- propmatch
}


plot(output,xlab="proportion of sites randomised",ylab="proportion of successful blast")
               


