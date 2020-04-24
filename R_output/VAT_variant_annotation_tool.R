require(VariantAnnotation)
require(ggplot2)
require(RJSONIO)
require(jsonlite)
require(httr)
require(plyr)
#make sure you have .vcf file in the same directory as working directory

VCF_file <- read.table(pipe("grep -v '^##' Challenge_data.vcf | sed s/^#//"),stringsAsFactors=F,header=T,sep="\t")

tempus_VCF_analysis <- VAT(VCF_file)

#First run VAT function and require packages

VAT <- function(VCF_file){
  
  VCF_matrix <- VCF_file

  VCF_annotated_file <- data.frame()
  for (i in 1:nrow(VCF_matrix)){
    
    #extracting informaion from Tempus Challenge File
    
    info <- as.character(VCF_matrix[i,]$INFO)
    info_collect <- sapply(info, function(x) strsplit(x, ";")[[1]], USE.NAMES=FALSE)
    info_fragment <- as.data.frame(t(sapply(info_collect, function(x) strsplit(x, "=")[[1]], USE.NAMES=FALSE)))
    info_VA <- setNames(as.character(info_fragment$V2), info_fragment$V1)
    
    #necessary info
    
    Reference_seq <- VCF_matrix[i,]$CHROM                     #the name of the sequence (i.e., chromosome)
    Position <- VCF_matrix[i,]$POS                            # 1-base position of the variation on the given sequence
    Reference_base <- VCF_matrix[i,]$REF                      # the reference base
    Alternative_base <- VCF_matrix[i,]$ALT                    #the list of alternative alleles at given position
    Reference_ID <- paste(Reference_seq,Position,
                          Reference_base,Alternative_base, 
                          sep = "-")                          #Reference ID for ExAC
    Variant_type <- info_VA[["TYPE"]]                         # variant type
    Read_depth <- as.numeric(info_VA[["DP"]])                 # Total read depth at the locus
    Allele_frequency <- as.numeric(info_VA[["AF"]])           #allele frequency for each ALT
    Allele_count <- as.numeric(info_VA[["AC"]])               # allele count in genotypes
    Number_of_allele <- as.numeric(info_VA[["AN"]])           # total number of alleles in genotypes
    Number_of_samples <- as.numeric(info_VA[["NS"]])          #number of samples in the data
    Ref_allele_quality <- as.numeric(info_VA[["QR"]])	        #Reference allele quality sum in phred
    Alt_allele_quality <- as.numeric(info_VA[["QA"]])	        #Alternate allele quality sum in phred
    Ref_forward_strand <- as.numeric(info_VA[["SRF"]])	      #Number of reference observations on the forward strand
    Ref_reverse_strand <- as.numeric(info_VA[["SRR"]])	      #Number of reference observations on the reverse strand
    Alt_forward_strand <- as.numeric(info_VA[["SAF"]])	      #Number of alternate observations on the forward strand
    Alt_reverse_strand <- as.numeric(info_VA[["SAR"]])	      #Number of alternate observations on the reverse strand
    Allele_balance <- as.numeric(info_VA[["AB"]])	            #Allele balance at heterozygous sites: a number between 0 and 1 representing the ratio of reads showing the reference allele to all reads
    Mean_mapping_quality_aa <- as.numeric(info_VA[["MQM"]])	  #Mean mapping quality of observed alternate alleles
    Mean_mapping_quality_ra <- as.numeric(info_VA[["MQMR"]])	#Mean mapping quality of observed reference alleles
    Proportion_of_aa <- as.numeric(info_VA[["PAIRED"]])	      #Proportion of observed alternate alleles which are supported by properly paired read fragments
    Proportion_of_ra <- as.numeric(info_VA[["PAIREDR"]])	    #Proportion of observed reference alleles which are supported by properly paired read fragments
    Ref_allele_obs_count <- as.numeric(info_VA[["RO"]])	      #Reference allele observation count
    Quality_ref_obs <- as.numeric(info_VA[["QR"]])	          #Sum of quality of the reference observations
    Alt_allele_obs_count <- as.numeric(info_VA[["AO"]])	      #Alternate allele observation count
    Quality_alt_obs <- as.numeric(info_VA[["QA"]])            #Sum of quality of the alternate observations
    
    #combining all the info
    
    Extracted_information <- cbind(Reference_ID,Reference_seq, Position,Reference_base, Alternative_base,
    Variant_type,Read_depth, Allele_frequency,Allele_count, Number_of_allele,
    Number_of_samples,Ref_allele_quality,Alt_allele_quality, Ref_forward_strand,
    Ref_reverse_strand,Alt_forward_strand,Alt_reverse_strand,  Allele_balance,
    Mean_mapping_quality_aa, Mean_mapping_quality_ra,Proportion_of_aa,
    Proportion_of_ra, Ref_allele_obs_count, Quality_ref_obs,Alt_allele_obs_count,
    Quality_alt_obs)
    
    #VCF annoted file to use further
    VCF_annotated_file <- rbind(VCF_annotated_file, Extracted_information)
    
  }
  

  
  # Connecting with Broad Institute ExAC Project API
  url = "http://exac.hms.harvard.edu/rest/bulk/variant"
  file_need <- toJSON(as.character(VCF_annotated_file$Reference_ID))
  
  ExAC_API_Content <- content(POST(url, body= file_need, encode = "json"))

  
  #VCF annotated file with adding neccessary info
  ExAC_Variant_Information <- vector()
  for (j in 1:nrow(VCF_annotated_file))
  {
    ID <- as.character(VCF_annotated_file$Reference_ID[j])
    ExAC_file <- ExAC_API_Content[[ID]]
    
    ExACinfo <- ExAC_API_Content[[as.character(VCF_annotated_file$Reference_ID[j])]]
    
    #gathering all necessary info
    
    Reference_ID <- as.character(VCF_annotated_file$Reference_ID[j])
    Gene_ID_ExAC <- paste("",as.character( ifelse(ExAC_file$variant$genes %in% c(""," ","NA"), NA, ExAC_file$variant$genes)),"",collapse=", ",sep="")
    #gene_symbol <-paste("",as.character(ifelse((ExAC_file$variant$vep_annotations[[1]]$SYMBOL) %in% c(""," ","NA"), NA, (ExAC_file$variant$vep_annotations[[1]]$SYMBOL))),"",collapse=", ",sep="")
    Transcript_ID_AC <-  paste("",as.character((as.character(ifelse(ExAC_file$variant$transcripts %in% c("NULL"," ","NA"), NA, ExAC_file$variant$transcripts)))),"",collapse=", ",sep="")
    Allele_freq_ExAC <-paste("",as.character(ifelse(ExAC_file$variant$allele_freq %in% c("NULL"," ","NA"), NA, ExAC_file$variant$allele_freq)),"",collapse=", ",sep="")
    Consequence_ExAC <- paste("",as.character(ifelse(names(ExAC_file$consequence) %in% c("NULL"," ","NA"), NA, names(ExAC_file$consequence))),"",collapse=", ",sep="")
    Reference_SNP_ID_ExAC <- paste("",as.character( ifelse(ExAC_file$variant$rsid %in% c("NULL"," ","NA"), NA, ExAC_file$variant$rsid)),"",collapse=", ",sep="")
    Population_info_ExAC <-paste("",as.character( ifelse(names(ExAC_file$variant$pop_ans) %in% c(""," ","NA"), NA, names(ExAC_file$variant$pop_ans))),"",collapse=", ",sep="")
   
     #combining all neccessary info
    
    EXAC_Info_combined <- cbind(Reference_ID,Gene_ID_ExAC,Transcript_ID_AC,Allele_freq_ExAC,
                                Allele_freq_ExAC, Consequence_ExAC,Reference_SNP_ID_ExAC,
                                Population_info_ExAC)
    
    ExAC_Variant_Information <- as.data.frame(rbind(ExAC_Variant_Information, EXAC_Info_combined))
  }
  
  #merging given vcf file info with ExAC info
  
  Annotated_VCF_file <- join(VCF_annotated_file, ExAC_Variant_Information,type = "inner")
  
  #save annotated VCF file
  save(Annotated_VCF_file, file="Annotated_VCF_file.rda") #save as R object
  write.table(Annotated_VCF_file, "Tempus_output_Annotated_VCF_file.txt", sep = "\t")
  
  
  return(Annotated_VCF_file)
  
  
  }

#Read Depth plot

depth_plot <- ggplot(tempus_VCF_analysis, aes(x=as.numeric(tempus_VCF_analysis$Read_depth))) + 
  geom_histogram(fill="gray54", alpha=0.5, bins = 100) +
  theme_bw() +xlab("Read Depth")+  ggtitle("Read Depth (DP)") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
depth_plot
pdf("read_depth_plot.pdf")
depth_plot
dev.off()


#Variant types 
variant_type_plot <- ggplot(tempus_VCF_analysis, aes(x = "", fill = factor(tempus_VCF_analysis$Variant_type))) + 
  geom_bar(width = 1) +
  coord_polar(theta = "y", start=0)+
  theme_bw()+ ggtitle("Variant Types") +
  labs(fill="Types of variants")+
  theme_void()
variant_type_plot
pdf("variant_type_plot.pdf")
variant_type_plot
dev.off()
