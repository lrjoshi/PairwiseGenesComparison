\name{pairwise_genes_comparsion}
\alias{pairwise_genes_comparsion}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
This package can be used to compare indiviual genes of the two genomes. This package makes use of coding sequences format that can be found in Genbank. Right now, the package allows pairwise comparison of the genes of two genomes. You find detail manual in the github wiki section. You will need FASTA Nuceotide sequence files of the two genomes. These can be obtained from the Genbank. If you click on Send to button on Genbank, then you can select "Coding Sequences" and then FASTA Nucleotide or FASTA Protein sequences. Then you can run pairwise_genes_comparison() funtion in R. You have to provide path for two files and then you can run pairwise_genes_comaprison function on those two files. It will create a nice summary table output. It will also save the table as a csv file called "summary_table,csv" in the current folder. It uses DECIPHER package for pairwise alignment and to generate similarity matrix.
}
\usage{
pairwise_genes_comparsion(filepath1, filepath2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  filepath1= path to the first genome file

  filepath2= path to the second genome file
}
\details{
Dependency : DECIPHER, TIDYVERSE, FORMATTABLE
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Lok Raj Joshi
}
\note{

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

#First create paths for the two files that you want to compare
filename1 <- paste("D:/Rpackages/gene1_sequence_format.txt")
filename2 <- paste("D:/Rpackages/gene2_sequence_format.txt")


#Then run pairwise_genes_comparison function using two filepaths
pairwise_genes_comparsion(filename1, filename2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
