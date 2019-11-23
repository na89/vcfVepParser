require(magrittr)
args = commandArgs(trailingOnly = TRUE)
inputFile <- args[1]
outputFile <- args[2]

##annotation rating
annotationRanking <-
  c(
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "transcript_amplification",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",
    "splice_region_variant",
    "incomplete_terminal_codon_variant",
    "start_retained_variant",
    "stop_retained_variant",
    "synonymous_variant",
    "coding_sequence_variant",
    "mature_miRNA_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
    "intron_variant",
    "NMD_transcript_variant",
    "non_coding_transcript_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "TFBS_ablation",
    "TFBS_amplification",
    "TF_binding_site_variant",
    "regulatory_region_ablation",
    "regulatory_region_amplification",
    "feature_elongation",
    "regulatory_region_variant",
    "feature_truncation",
    "intergenic_variant"
  )

incon <- gzcon(file(inputFile, open="rb"))
outConAll <- file(paste(outputFile, "_all", sep = ""), "w")
outConSevere <- file(paste(outputFile, "_severe", sep = ""), "w")

while ( TRUE ) {
  line = readLines(incon, n = 1)
  if ( length(line) == 0 ) {
    break
  }
  if ( grepl(pattern = "##", x = line) ){
    next
  }
  originalLine <- line
  line <- strsplit(x = line, split = "\t") %>%
    unlist %>%
    strsplit(x = .[8], split = ";")
  
  csqPosition <- lapply(line, function(x){
    grepl(pattern = "CSQ=", x = x)
  }) %>%
    unlist
  csqPosition <- which(csqPosition)
  line <- line[[csqPosition]] %>%
    unlist
  line <- strsplit(x = line, split = ",") %>%
    unlist
  effect <- lapply(line, function(x){
    strsplit(x, split = "\\|") %>%
      unlist %>%
      .[2]
  }) %>%
    unlist
  
  gene <- lapply(line, function(x){
    strsplit(x, split = "\\|") %>%
      unlist %>%
      .[4]
  }) %>%
    unlist
  
  outputPrefix <- strsplit(x = originalLine, split = "\t") %>%
    unlist %>%
    .[c(1,2,4,5)] %>%
    paste(., collapse = "\t")
  
  if(length(effect) != length(gene)){
    message("Error at line ", originalLine)
    stop()
  }
  for (i in 1:length(effect)){
    writeLines(paste(outputPrefix, effect[i], gene[i], sep = "\t"), 
               con = outConAll)
  }
  severeImpact <- min(which(annotationRanking %in%
                          effect))
  severeImpact <- annotationRanking[severeImpact]
  severeGene <- gene[which(effect == severeImpact)]
  
  for(i in 1:length(severeGene)){
    writeLines(paste(outputPrefix, severeImpact, severeGene[i], sep = "\t"),
               con = outConSevere)
  }
  
  
}

close(incon)
close(outConSevere)
close(outConAll)


