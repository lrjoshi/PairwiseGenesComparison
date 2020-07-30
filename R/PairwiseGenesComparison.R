#required packages
library (tidyverse)
library (readr)
library (seqinr)
library (DECIPHER)
library (formattable)


#Provide path to the filename
#filename1 <- paste("D:/Courses/NGS Tutorials/R_Practice/gff3/gene_sequence_format.txt")
#filename2 <- paste("D:/Courses/NGS Tutorials/R_Practice/gff3/gene_sequence_format.txt")

#filename1 <- paste ("SARS-Cov-2-Italy.txt")
#filename2<- paste("SARS-CoV-2-Wuhan.txt")

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


  #extracting gene name ,start and end position


  #Change identifier for file 1 here "gene" can be replaced with locus tag, protein, protein ID.
  table1 <- table1 %>%
    mutate(name=  str_replace_all(name, "^.*?gene=(.*?)\\].*?\\[location.*?(\\d+)\\.\\.(\\d+).*?$",
                                  "\\1___\\2___\\3")) %>%
    separate(name,
             sep="___",
             into = c("gene", "start", "end"))


  #Change identifier for file 2 here "gene" can be replaced with locus tag, protein, protein ID.

  table2 <- table2 %>%
    mutate(name=  str_replace_all(name, "^.*?gene=(.*?)\\].*?\\[location.*?(\\d+)\\.\\.(\\d+).*?$",
                                  "\\1___\\2___\\3")) %>%
    separate(name,
             sep="___",
             into = c("gene", "start", "end"))






  #### ROW WISE ORGANIZATION OF GENES ####
  #create pairwise files
  #it cretes CSV files in the pairs directory
  #one file for each gene
  #rowwise organization of genes
  for (i in 1:nrow(table1)){
    table3 <- rbind(table1[i,],table2[i,])
    table3 <- as.data.frame(table3)
    table3<-select(table3,2:5)
    filename <- paste0("genes",i,".csv")
    write.csv(table3,filename,row.names = F)
  }


  #### CONVERT TABLE INTO FASTA ####
  #this will convert table into fasta file
  #one fasta file for a pair of gene
  #apply tabular to fasta funtion to all the dataframes



  for (j in 1:nrow(table1)){
    infilename = paste0("genes",j,".csv")
    file <- read.csv (file=infilename, header=TRUE)

    file = as.data.frame(file)


    #file$gene=c(name_vector1[j],name_vector2[j])

    #delete if any existing file

    # unlink(c("dna_fasta.fasta"), force=TRUE)

    #give output filename

    outfilename = paste0("genes",j,".fasta")
    sink(outfilename)

    for (k in 1:nrow(file)){
      name = paste0(">",file[k,1])
      sequence = paste(file[k,4])
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


  #format similarity to two decimal places
  percent_similarity=(1-similarity)*100
  percent_similarity= format(round(percent_similarity, 2), nsmall = 2)

  #create the table of the genes and the similarity score between two genes

  similarity_table <- data.frame (Gene_Name1=table1$gene,
                                  Start_Sequence1=table1$start,
                                  Stop_Sequence1=table1$end,
                                  Gene_Length1= abs(as.numeric(table1$start)-as.numeric(table1$end)),
                                  Gene_Name2 = table2$gene,
                                  Start_Sequence2=table2$start,
                                  Stop_Sequence2=table2$end,
                                  Gene_Length2= abs(as.numeric(table2$start)-as.numeric(table2$end)),
                                  Percentage_Similarity=percent_similarity)





  write.csv(similarity_table,"summary_table.csv")



  #### creating nice table ####
  df_format <- read_csv ("summary_table.csv")

  nice_table <- formattable (df_format,align=c("c","c","c","c","c","c","c","c","c","l"),
                             list(
                               Gene_Name1= formatter("span",style = ~ style(color = "black",font.weight = "bold")),
                               Gene_Name2= formatter("span",style = ~ style(color = "black",font.weight = "bold")),
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


#filename1 <- paste ("SARS-Cov-2-Italy.txt")
#filename2<- paste("SARS-CoV-2-Wuhan.txt")
#pairwise_genes_comparison(filename1,filename2)
