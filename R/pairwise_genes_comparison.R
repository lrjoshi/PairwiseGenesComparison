#required packages
library (tidyverse)
library (readr)
library (seqinr)
library (DECIPHER)
library (formattable)


#Provide path to the filename
filename1 <- paste("D:/Courses/NGS Tutorials/R_Practice/gff3/gene_sequence_format.txt")
filename2 <- paste("D:/Courses/NGS Tutorials/R_Practice/gff3/gene_sequence_format.txt")



pairwise_genes_comparison <- function (filename1,filename2) {


file1 <- readLines(filename1)

file2 <- readLines(filename2)


#### Parsing Files and Covert to CSV ####
#read fasta file1



#find the genename location by grepping >

location <- which((str_sub(file1,1,1))==">")

#start an empty vector to collect name and sequence

name=c()
sequence =c()



#number of genes= number of loops
#extract name first
for ( i in 1:length(location)){
  name_line = location[i]
  name1 = file1[name_line]
  name=c(name,name1)
  #extract sequence between the names
  #the last sequence will be missed using this strategy
  #so, we are using if condition to extract last sequence
  start= location[i]+1
  end = location[i+1]-1
  if ( i < length (location)){

    end=end

  } else {

    end=length(file1)
  }

  lines = start:end
  sequence1= as.character(paste(file1[lines],collapse = ""))
  sequence =c(sequence,sequence1)
}

#now create table using name and sequence vector

data <- tibble(name=name,sequence=sequence)



#finally export the file
#before that remove preexisting file
unlink(c("dna_table1.csv"),force=TRUE)
write.csv(data,"dna_table1.csv")




#read fasta file2



#find the genename location by grepping >

location2 <- which((str_sub(file2,1,1))==">")

#start an empty vector to collect name and sequence

name2=c()
sequence2 =c()



#number of genes= number of loops
#extract name first
for ( i in 1:length(location2)){
  name_line2 = location2[i]
  nam = file2[name_line2]
  name2=c(name2,nam)
  #extract sequence between the names
  #the last sequence will be missed using this strategy
  #so, we are using if condition to extract last sequence
  start2= location2[i]+1
  end2 = location2[i+1]-1
  if ( i < length (location2)){

    end2=end2

  } else {

    end2=length(file2)
  }

  lines2 = start2:end2
  sequen= as.character(paste(file2[lines2],collapse = ""))
  sequence2 =c(sequence2,sequen)
}

#now create table using name and sequence vector

data <- tibble(name=name2,sequence=sequence2)



#finally export the file
#before that remove preexisting file
unlink(c("dna_table2.csv"),force=TRUE)
write.csv(data,"dna_table2.csv")









#data <- paste("gene_sequence_format.txt")


table1 <- read.csv("dna_table1.csv",header=T)
table2 <- read.csv("dna_table2.csv",header=T)


#### Formatting the table  ####
#extracing information from the name column
names1 <- table1$name
names2 <- table2$name

string_list1 <-  str_extract_all(names1, boundary("word"))
string_list2 <-  str_extract_all(names2, boundary("word"))



name_vector1 = c()
name_vector2 = c()
start_position1 =c()
start_position2 =c()
stop_position1 =c()
stop_position2 =c()
gene_length1=c()
gene_length2=c()
for (z in 1:nrow(table1)) {
  gene_name1 <- string_list1[[c(z,4)]]
  start_pos1 <- as.numeric(string_list1[[c(z,11)]])
  stop_pos1 <-as.numeric( string_list1[[c(z,12)]])


  gene_name2 <- string_list2[[c(z,4)]]
  start_pos2 <- as.numeric(string_list2[[c(z,11)]])
  stop_pos2 <- as.numeric(string_list2[[c(z,12)]])


  gene_length2 =append(gene_length2,(stop_pos2-start_pos2+1))
  gene_length1 =append(gene_length1,(stop_pos1-start_pos1+1))

  name_vector1 <- append(name_vector1,gene_name1)
  start_position1 <- append (start_position1,start_pos1)
  stop_position1 <- append (stop_position1,stop_pos1)

  name_vector2 <- append(name_vector2,gene_name2)
  start_position2 <- append (start_position2,start_pos2)
  stop_position2 <- append (stop_position2,stop_pos2)
}


#check if names are similar
#otherwise terminate the funtion
#if gene names are not identical the analysis will not be correct


validation <- identical(name_vector1,name_vector2)


if (validation==T) {
  print("Validation complete. Continuing analysis.")
} else if (validation== F) {
  stop("Stopping analysis. The name, number or order of the genes in the selected files are not same. Check package docmentation for details.")
}




#### ROW WISE ORGANIZATION OF GENES ####
#create pairwise files
#it cretes CSV files in the pairs directory
#one file for each gene
#rowwise organization of genes
for (i in 1:nrow(table1)){
  table3 <- rbind(table1[i,],table2[i,])
  table3 <- as.data.frame(table3)
  table3<-select(table3,2:3)
  filename <- paste0("genes",i,".csv")
  write.csv(table3,filename,row.names = F)
}


#### CONVERT TABLE INTO FASTA ####
#this will convert table into fasta file
#one fasta file for a pair of gene
#apply tabular to fasta funtion to all the dataframes



for (j in 1:11){
  infilename = paste0("genes",j,".csv")
  file <- read.csv (file=infilename, header=TRUE)

  file = as.data.frame(file)
  #delete if any existing file

  # unlink(c("dna_fasta.fasta"), force=TRUE)

  #give output filename
  outfilename = paste0("genes",j,".fasta")
  sink(outfilename)

  for (k in 1:nrow(file)){
    name = paste0(">",file[k,1])
    sequence = paste(file[k,2])
    cat(name,sep="\n")
    cat(sequence,sep="\n")
  }

  #this is sink all the console output to the file
  sink()

}




#### CREATING ALIGNEMENT WITH DECIPHER ####
###now creating alignment of the files
#get the similarity of the pairwise alignment
#print the alignement files


similarity =c()
for (l in 1:nrow(table1)){
  fasta_filename <- paste0("genes",l,".fasta")
  seqs<-readDNAStringSet(fasta_filename)
  aligned <- AlignSeqs(seqs)
  #print(aligned)
  d<- DistanceMatrix(aligned)
  distance <- d[2]
  similarity <- append(distance,similarity)
  align_filenames <- paste0("aligned",l,".fasta")
  writeXStringSet(aligned,file=align_filenames)
}



#create the table of the genes and the similarity score between two genes

similarity_table <- data.frame (Gene_Name=name_vector1,
                                Start_Sequence1=start_position1,
                                Stop_Sequence1=stop_position1,
                                Gene_Length1= gene_length1,
                                Start_Sequence2=start_position2,
                                Stop_Sequence2=stop_position2,
                                Gene_Length2= gene_length2,
                                Percentage_Similarity=((1-similarity)*100))





write.csv(similarity_table,"summary_table.csv")

#### creating nice table ####
df_format <- read_csv ("summary_table.csv")

nice_table <- formattable (df_format,align=c("c","c","c","c","c","c","c","c","l"),
             list(
               Gene_Name= formatter("span",style = ~ style(color = "black",font.weight = "bold")),
               Start_Sequence1=color_tile ("#E4F9E3","#E4F9E3"),
               Stop_Sequence1=color_tile ("#E4F9E3","#E4F9E3"),
               Gene_Length1=color_tile ("#E4F9E3","#E4F9E3"),
               Start_Sequence2=color_tile ("#D1E6FA","#D1E6FA"),
               Stop_Sequence2=color_tile ("#D1E6FA","#D1E6FA"),
               Gene_Length2=color_tile ("#D1E6FA","#D1E6FA"),

               Percentage_Similarity =color_bar ("#FA614B")
             ))

print (nice_table)
########Final Message#################

print ("Check summary_table.csv file in the current directory.")
print ("The suffix 1 and 2 denotes the file1 and file2 respectively.")
print ("Thank you. End of analysis.")

}
head(similarity_table)

pairwise_genes_comparison(filename1,filename2)
