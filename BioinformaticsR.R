#Download the data
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part1_files/mydata.tsv",
              destfile = "mydata.tsv")
              
#Read into R making sure the 1st column is set as the row names and headers are recognised
x <- read.table("mydata.tsv", header = TRUE , stringsAsFactors = FALSE )

#Lets check if the data has been imported properly
head(x)
str(x)
