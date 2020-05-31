library("seqinr")
library("R.utils")
library("rBLAST")
library("ape")
library("ORFik")
library("Biostrings")
#To download the e.coli file
download.file("ftp://ftp.ensemblgenomes.org/pub/bacteria/release-42/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/cds/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz",
              destfile = "e.coli.fa.gz")
R.utils::gunzip("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz", overwrite= TRUE)
unlink("e.coli.fa")
gunzip("e.coli.fa.gz")

makeblastdb("e.coli.fa",dbtype="nucl","-parse_seqids")
#To download sample file
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/sample.fa",
              destfile = ("sample.fa"))
x <- read.fasta("sample.fa")
myseq <- x[[2]]
myseq
str(myseq)
seqinr::getLength(myseq)
length(myseq)
seqinr::GC(myseq)

download.file("https://github.com/markziemann/SLE712_files/blob/master/bioinfo_asst3_part2_files/mutblast_functions.R",
              destfile = "mutblast.R")





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
tophit <- db[which(names(db) %in% hits[1])]
tophit[1:3]
res


mutator <- function(myseq,nmut) {
  myseq_mod <- myseq
  mypos<-sample(seq_along(myseq),nmut)
  myseq_mod[mypos] <- sample(c("a","c","g","t"),length(mypos),replace = TRUE)
  return(myseq_mod)
}






