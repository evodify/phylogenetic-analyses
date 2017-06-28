# Reconstructs NJ phylogenetic trees from many fasta files
# usage:
# Rscript NJ_slidingWindow.R *.fasta

library(ape)

file.create('trees.nwk')  # creat empty file
file.create('trees.names')  # creat empty file

ff <- commandArgs(trailingOnly = T) # a list of files to process

for (file in ff) {
  snps <- read.dna(file, format = "fasta")
  dd <- dist.dna(snps, model ='raw', pairwise.deletion = TRUE)
  tryCatch({
    NJtree <- nj(dd)
    NJtreeR <- root(NJtree, 'Neslia')
    NJtreeR$edge.length <- NULL # remove branch length
    write.tree(NJtreeR, 'trees.nwk', append=T)
    write(file, 'trees.names', append=T)
    message(paste(file, 'processed', sep=' '))
    },
    error=function(e){message(paste(file, 'skipped', sep=' '))})
}

