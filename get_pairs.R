library(biomaRt)
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org")
ensembl_human <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
gene_info <- getBM(attributes = c("ensembl_gene_id",     
                                  "strand",              
                                  "chromosome_name",    
                                  "gene_biotype",         
                                  "start_position",      
                                  "end_position",        
                                  "transcription_start_site"),  
                   mart = ensembl_human)


locate_BIP <- function(gene_info){
  colnames(gene_info) <- c("Gene_ID", "Strand", "Chromosome", "Biotype",
                           "Gene_Start", "Gene_End", "TSS")
  
  gene_info <- gene_info[gene_info$Chromosome %in% c(1:22, "X"),]
  gene_info <- gene_info[order(gene_info$TSS),]
  gene_info <- gene_info[order(gene_info$Chromosome),]
  gene_info <- gene_info[!duplicated(gene_info$Gene_ID),]
  
  ##Add columns for distance and pair IDs
  gene_info$Distance <- NA
  gene_info$PairID <- NA
  
  ##assign input data frame length to len and count to 0
  len <- dim(gene_info)[1]
  count <- 0
  
  ##initialise empty dataframe in which to append putative gene pairs
  putative_bidirectional <- data.frame()
  
  ##loop through gene_info and compare the TSS to the TSS of the next gene in
  ##the data frame.
  for (row in 1:(len-1)){
    
    ##assign variables
    A <- gene_info$TSS[row]
    B <- gene_info$TSS[row+1]
    strandA <- gene_info$Strand[row]
    strandB <- gene_info$Strand[row+1]
    chromeA <- gene_info$Chromosome[row]
    chromeB <- gene_info$Chromosome[row+1]
    BiotypeA <- gene_info$Biotype[row]
    BiotypeB <- gene_info$Biotype[row+1]
    dist <- abs(A-B)
    
    ##update distance column
    gene_info$Distance[row] <- dist
    
    ##criteria:
    ##less than 1000bp apart
    ##on opposite strands 
    ##further than 100bp apart
    ##on the same chromosome
    ##the first in the pair is on the -1 strand 
    ##(and thus excludes overlapping genes)
    ## & BiotypeA != BiotypeB & BiotypeA %in% c('lncRNA', 'protein_coding') & BiotypeB %in% c('lncRNA', 'protein_coding') ?
    if (dist < 1000 & strandA != strandB & dist > 100 & chromeA == chromeB &
        strandA == -1){
      
      ##increase count
      count <- count + 1
      
      ##append bidirectional pair to dataframe
      putative_bidirectional <- rbind(putative_bidirectional, gene_info[row,])
      putative_bidirectional <- rbind(putative_bidirectional, gene_info[row+1,])
      
      ##assign the index of the last entry in the bidirectional gene list, which
      ##was just added
      index <- as.integer(dim(putative_bidirectional)[1])
      
      ##use that index to update the pair ID
      padded_count <- formatC(count, width = 4, format = "d", flag = "0")
      putative_bidirectional$PairID[index-1] <- paste("Pair", as.character(padded_count), sep="_")
      putative_bidirectional$PairID[index] <- paste("Pair", as.character(padded_count), sep="_")
    }
  }
  putative_bidirectional$Distance <- putative_bidirectional$Distance - 1
  return(putative_bidirectional)
}

putative_bidirectional <- locate_BIP(gene_info)
