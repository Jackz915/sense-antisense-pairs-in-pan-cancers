
# Run jaspar
# ./jaspar_tfbs_extraction/bin/extract_TFBSs_JASPAR.sh -i putative_bidirectional.bed -b JASPAR2024_hg38.bb -s 600 -p 48 > output
# less -S output  | awk -F '\t' '{print $1"\t"$2"\t"$3"\t"$7}' | less  -S | sort | uniq  > TF_binding_sites

TF_binding_sites <- read.table('GDCdata/TF_binding_sites',header = F,sep = "\t")
colnames(TF_binding_sites) <- c('Chromosome','Start', 'End','TF')
TF_binding_sites <- TF_binding_sites %>%
  mutate(TF = strsplit(as.character(TF), "[::|-]")) %>%  
  unnest(TF)  

TF_binding_sites <- TF_binding_sites[TF_binding_sites$TF != '',]


putative_bidirectional <- putative_bidirectional %>%
  rowwise() %>%
  mutate(
    putative_TF = paste(
      unique(TF_binding_sites$TF[
        TF_binding_sites$Chromosome == Chromosome & 
          ((TF_binding_sites$Start >= TSS_sense & TF_binding_sites$Start <= TSS_antisense) |
             (TF_binding_sites$End >= TSS_sense & TF_binding_sites$End <= TSS_antisense) |
             (TF_binding_sites$Start <= TSS_sense & TF_binding_sites$End >= TSS_antisense))]),
      collapse = ","
    )
  ) %>%
  ungroup()
