# Loads all required libraries and verifies that they are present.

libload <- c("assertthat", "stringr", "tools", "IRanges", "GenomicAlignments", "ggrepel", "dplyr", "Rsamtools", "Rsubread", "DESeq2", "knitr", "rmarkdown", "testthat", "baerhunter", "GenomicRanges", "rtracklayer")
lapply(libload, require, character.only = TRUE)

# Sets the working directory of the project
setwd("/home/princeofparthia/Documents/mtb_rna")

# The original feature_file_editor function from Baerhunter
feature_file_editor(bam_directory="PRJNA480455/test/",  
                    original_annotation_file="original_annotation.gff3",
                    annot_file_dir = "PRJNA480455/gff3/",
                    output_file="PRJNA480455/gff3/test.gff3", 
                    original_sRNA_annotation="ncRNA", 
                    low_coverage_cutoff=5, 
                    min_sRNA_length=40, 
                    high_coverage_cutoff=10, 
                    min_UTR_length=50, 
                    paired_end_data=TRUE, 
                    strandedness="stranded")

# The original count_features function from Baerhunter.
count_features(bam_dir = "PRJNA480455/test/",
               annotation_dir= "PRJNA480455/gff3/", 
               annotation_file = "test.gff3", 
               output_dir = "PRJNA480455/test/",
               chromosome_alias_file = "PRJNA480455/aliasfile.txt",
               strandedness = "stranded", 
               is_paired_end= TRUE)

# The wrapper function coverage_cutoff_iterator around the Baerhunter functions with customisation and the ability to step through coverage cutoffs.

#' @param bam_directory The directory containing BAM files.
#' @param original_annotation_file The name of the referenc GFF3 genome annotation file.
#' @param annot_file_dir The directory containing the reference GFF3 annotation file.
#' @param output_file_name A variable containing the name of a GFF3 output file. Changes the output file name for each step to allow saving of multiple steps from the same run. Can be a customised string.
#' @param low_coverage_start An integer indicating the first low coverage threshold value to test. Low coverage is the minimum value of coverage for this function to start testing for putative ncRNA
#' @param low_coverage_end An integer indicating the last low coverage threshold value to test.
#' @param low_coverage_search_factor An integer indicating by what amount the low coverage threshold will increase per iteration.
#' @param high_coverage_start An integer indicating the first high coverage threshold value to test. High coverageis the minimum value the coverage needs to maintain in order to be counted as putative ncRNA
#' @param high_coverage_end An integer indicating the last high coverage threshold value to test.
#' @param high_coverage_search_factor An integer indicating by what amount the high coverage threshold will increase per iteration.
#' @param dataset_name A variable containing the name of the dataframe file produced with all ncRNA counts. Can be a customised string.
#' @param output_dir The directory to output all GFF3 files to. This needs to be created before starting the function.
#' @param count_file_dir The directory to output all Counts and Counts_summary files to.
#' @param moving_frame A boolean variable to determine if the low coverage and high coverage cutoffs move in tandem (e.g 1000 - 1200, 1200 - 1400, etc.)
#' @param strandedness A string indicating the strandedness of the dataset. Can be "stranded", "reversely stranded" and "unstranded".

coverage_cutoff_iterator <- function(bam_directory = ".", 
                                   original_annotation_file = ".gff3", 
                                   annot_file_dir = "", 
                                   output_file_name = standard,
                                   low_coverage_start = 10, 
                                   low_coverage_end = 50,
                                   low_coverage_search_factor = 10,
                                   high_coverage_start = 100, 
                                   high_coverage_end = 150,
                                   high_coverage_search_factor = 10,
                                   dataset_name = standard,
                                   output_dir = "", 
                                   chromosome_alias_file = "", 
                                   count_file_dir = "",
                                   moving_frame = TRUE,
                                   strandedness = "stranded") {
  
  # Sets number of "steps" that need to occur in each pass (hc), and the number of passes (lc)
  lc = (low_coverage_end - low_coverage_start)/low_coverage_search_factor
  hc = (high_coverage_end - high_coverage_start)/high_coverage_search_factor
  # Creates list variables
  directory <- c()
  low_coverage_cutoff <- c()
  high_coverage_cutoff <- c()
  sRNA <- c()
  UTR <- c()
  # Some simple checks to ensure coverage_ends are larger than coverage_starts.
  if (lc<0) {
    message("low-coverage cutoff start larger than low-coverage cutoff end!")
    break
  }
  if (hc<0){
    message("high-coverage cutoff start larger than high-coverage cutoff end!")
    break
  }
  # Calculates total number of iterations to be printed in message() (this will be wrong if moving_frame = TRUE or if low_coverage_end is higher than high_coverage_start, as many steps will be skipped.)
  k <- (lc+1)*(hc+1)
  for (i in 0:lc) {
    # Calculates the low and high_coverage cutoff values to be used in each individual step, for reporting to the user and saving GFF3 files.
    x = low_coverage_start + i*low_coverage_search_factor
    for (j in 0:hc) {
      y = high_coverage_start + j*high_coverage_search_factor
      # Messages the user  for all instances where the high coverage cutoff is equal to or below the low coverage cutoff, and skips that iteration.
      if (y<=x) {
        message("Skipping ", (i+1), " low coverage tests, ", (j+1), " high coverage tests, out of ", k, " total tests, or ", x, ", ", y, ": low coverage cutoff is larger than or equal to high coverage cutoff.")
        next
      }
      # Causes moving frame to skip all iterations where the high coverage cutoff is further away from the low coverage cutoff than one step of the high coverage cutoff. This allows both cutoffs to move in tandem.
      if (moving_frame == TRUE){
         if (y-high_coverage_search_factor>x) {
          message("Skipping ", (i+1), " low coverage tests, ", (j+1), " high coverage tests, out of ", k, " total tests, or ", x, ", ", y, ": (moving frame).")
          next
        } 
      }
      # Allows customisation of the output_file_name variable.
      if (output_file_name == standard) {
      output_file= paste0(output_dir, "/Coverage_", x, "_", y, ".gff3")
      feature_file_editor(bam_directory= bam_directory,  
                          original_annotation_file=original_annotation_file,
                          annot_file_dir = annot_file_dir,
                          output_file=output_file, 
                          original_sRNA_annotation="ncRNA", 
                          low_coverage_cutoff= x, 
                          min_sRNA_length=30, 
                          high_coverage_cutoff=y, 
                          min_UTR_length=50, 
                          paired_end_data=TRUE, 
                          strandedness=strandedness)
      }
      else {
        output_file= paste0(output_dir, "/Coverage_", x, "_", y, "_", output_file_name, ".gff3")
        feature_file_editor(bam_directory= bam_directory,  
                            original_annotation_file=original_annotation_file,
                            annot_file_dir = annot_file_dir,
                            output_file=output_file, 
                            original_sRNA_annotation="ncRNA", 
                            low_coverage_cutoff= x, 
                            min_sRNA_length=30, 
                            high_coverage_cutoff=y, 
                            min_UTR_length=50, 
                            paired_end_data=TRUE, 
                            strandedness=strandedness)
      }
      # Allows customisation of the dataset_name variable.
      if (dataset_name == standard) {  
        count_features(bam_dir = bam_directory,
                   annotation_file = output_file, 
                   output_dir = "PRJNA480455/test/",
                   output_filename = paste0(x, "_", y, "_dataset"),
                   chromosome_alias_file = chromosome_alias_file,
                   strandedness = strandedness, 
                   is_paired_end= TRUE)
      }
      else {
        count_features(bam_dir = bam_directory,
                       annotation_file = output_file, 
                       output_dir = "PRJNA480455/test/",
                       output_filename = dataset_name,
                       chromosome_alias_file = chromosome_alias_file,
                       strandedness = strandedness,
                       is_paired_end= TRUE)
      }
  # Script modified from merging GFF3 file step in order to create a dataframe of coverage cutoff windows and the ncRNA elements Baerhunter predicts within them.    
  g1 <- import.gff3(output_file)
  
  srnas_gr1<-g1[elementMetadata(g1)[,"type"] %in% "putative_sRNA"]
  srnas_gr1<-srnas_gr1[which(width(srnas_gr1)<=1000),]
  utrs_gr1<-g1[elementMetadata(g1)[,"type"] %in% "putative_UTR"]
  utrs_gr1<-utrs_gr1[which(width(utrs_gr1)<=1000),]
  
  # Appends all data from this iteration of the function in each relavent list.  
  directory <- append(directory, bam_directory)
  low_coverage_cutoff <- append(low_coverage_cutoff, x)
  high_coverage_cutoff <- append(high_coverage_cutoff, y)
  sRNA <- append(sRNA, length(srnas_gr1))
  UTR <- append(UTR, length(utrs_gr1))
  # Messages a completetion step when that coverage cutoff window has finished.
  message((i+1), " low coverage tests, ", (j+1), " high coverage tests, out of ", k, " total tests done")
  }
  }
  # Writes the dataframe from the lists of each data type from all iterations. 
  df <- DataFrame(directory, low_coverage_cutoff, high_coverage_cutoff, sRNA, UTR)
  # Writes dataframe to file.
  write.csv(df, count_file_dir, row.names=FALSE)
  # Also returns the dataframe for further use in R / inspection by user.
  return(df)
}

# Example use case of coverage_cutoff_iterator
df <- coverage_cutoff_iterator(bam_directory = "PRJNA480455/glucose", 
                       original_annotation_file = "original_annotation.gff3",
                       annot_file_dir = "PRJNA480455/gff3/",
                       low_coverage_start = 1000,
                       low_coverage_end = 8000,
                       low_coverage_search_factor = 200,
                       high_coverage_start = 1200,
                       high_coverage_end = 8000,
                       high_coverage_search_factor = 200,
                       output_dir = "PRJNA480455/test/glucose_output",
                       chromosome_alias_file = "PRJNA480455/aliasfile.txt",
                       count_file_dir = "PRJNA480455/test/glucose_output")


### Script from J. Stiens' code to combine the GFF3 files produced above into a combined GFF3 file. Modified to work with this dataset.

# adapted from script comb_gff.R
# combining (filtered) gff3 files with predicted sRNA/UTR annotations to get one master list

# 1. Subset elements by ncRNAs (sRNA and UTR) 
# 2. make genomic ranges object for each treatment. filter for length
# 3. combine all ncRNAs together in one grange object
# 4. reduce to simplified set of ncRNA coordinates
# 5. map back to get names of ncRNAs? if combine sRNAs and UTRs, to identify which? 
# import each .gff and subset by ncRNAs

ref <- import.gff3("PRJNA480455/gff3/original_annotation.gff3")

g1 <- import.gff3("PRJNA480455/test/Coverage_6000_6200_all_glucose.gff3")
g2 <- import.gff3("PRJNA480455/test/Coverage_6000_6200_all_lactate.gff3")
g3 <- import.gff3("PRJNA480455/test/Coverage_6000_6200_all_pyruvate.gff3")

# add all sRNAs and utrs to common list (do these separately)

srnas_gr1<-g1[elementMetadata(g1)[,"type"] %in% "putative_sRNA"]

# filter for sRNAs/UTRs > 1000 nt

srnas_gr1<-srnas_gr1[which(width(srnas_gr1)<=1000),]
utrs_gr1<-g1[elementMetadata(g1)[,"type"] %in% "putative_UTR"]
utrs_gr1<-utrs_gr1[which(width(utrs_gr1)<=1000),]

srnas_gr2<-g2[elementMetadata(g2)[,"type"] %in% "putative_sRNA"]
srnas_gr2<-srnas_gr2[which(width(srnas_gr2)<=1000),]
utrs_gr2<-g2[elementMetadata(g2)[,"type"] %in% "putative_UTR"]
utrs_gr2<-utrs_gr2[which(width(utrs_gr2)<=1000),]

srnas_gr3<-g3[elementMetadata(g3)[,"type"] %in% "putative_sRNA"]
srnas_gr3<-srnas_gr3[which(width(srnas_gr3)<=1000),]
utrs_gr3<-g3[elementMetadata(g3)[,"type"] %in% "putative_UTR"]
utrs_gr3<-utrs_gr3[which(width(utrs_gr3)<=1000),]

total_srnas_gr<-c(srnas_gr1, srnas_gr2, srnas_gr3)
length(total_srnas_gr)

# reduce to align ranges and merge overlapping ranges, those with specified gap between 

red_srnas_gr<-reduce(total_srnas_gr, min.gapwidth=5L)
length(red_srnas_gr)

# do the same with utrs

total_utrs_gr<-c(utrs_gr1, utrs_gr2, utrs_gr3)
length(total_utrs_gr)

red_utrs_gr<-reduce(total_utrs_gr, min.gapwidth=5L)
length(red_utrs_gr) 

# merge overlapping srnas and utrs to find srnas that overlap utr ranges

##mergeByOverlaps(query, subject, ...)"

ov <- mergeByOverlaps(red_srnas_gr, red_utrs_gr)

# for each element in overlap, keep longer one

delete_utr <- GRanges()
delete_srna <- GRanges()

for (i in 1:nrow(ov)) {
  if (width(ov$red_srnas_gr)[i] > width(ov$red_utrs_gr)[i]) {
    delete_utr <- c(delete_utr, ov$red_utrs_gr[i])
  }
  else {
    delete_srna <- c(delete_srna, ov$red_srnas_gr[i])
  }
}

# remove ranges in red_srnas_gr that are in ov and smaller than utr

red_srnas_gr <- red_srnas_gr[!red_srnas_gr %in% delete_srna]
length(red_srnas_gr)

#remove ranges in red_utrs_gr that overlap srnas and are smaller

red_utrs_gr <- red_utrs_gr[!red_utrs_gr %in% delete_utr]
length(red_utrs_gr)

#save object

#saveRDS(red_srnas_gr, here("R_Data/reduced_srna_granges.RData"))
#saveRDS(red_utrs_gr, here("R_Data/reduced_utr_granges.RData"))

# read objects
#red_srnas_gr <-readRDS(here("R_Data/reduced_srna_granges.RData"))
#red_utrs_gr <- readRDS(here("R_Data/reduced_utr_granges.RData"))

# read in annotated sRNAs (type=ncRNA_gene) from gff

gff <- import("PRJNA480455/gff3/original_annotation.gff3")
gff_df <- read.delim("PRJNA480455/gff3/original_annotation.gff3", header=F, comment.char = "#")

#unique(gff_df$type)
#nrow(gff_df[gff_df$type=="biological_region",])

nc_gr <- gff[(elementMetadata(gff)[, "type"]=="ncRNA_gene")]

# filter out overlapping ranges with known annotated ncRNAs (utrs and srnas)

remove_srnas <- subsetByOverlaps(red_srnas_gr, nc_gr, ignore.strand=FALSE)

#9 ranges overlap annotated ncRNAs

remove_utrs  <- subsetByOverlaps(red_utrs_gr, nc_gr, ignore.strand=FALSE)
length(remove_utrs)

#13  (good since there are 30 ncRNA, covers 22 of them)

red_srnas_gr <- red_srnas_gr[!red_srnas_gr %in% remove_srnas]
red_utrs_gr  <- red_utrs_gr[!red_utrs_gr %in% remove_utrs]

#save object

#saveRDS(red_srnas_gr, here("R_Data/reduced_srna_granges.RData"))
#saveRDS(red_utrs_gr, here("R_Data/reduced_utr_granges.RData"))

# subset by strand

pos_srnas_gr<-subset(red_srnas_gr, strand=="+")
neg_srnas_gr<-subset(red_srnas_gr, strand=="-")

# this uses code from BH srna_calc and utr_calc to re-name the new elements 
# (needed to add 'as.integer' to remove whitespace before coordinates)

names(pos_srnas_gr) <- apply(as.data.frame(pos_srnas_gr), 1, function(x) paste("ID=putative_sRNA:p", as.integer(x[2]), "_", as.integer(x[3]), ";", sep = ""))
names(neg_srnas_gr) <- apply(as.data.frame(neg_srnas_gr), 1, function(x) paste("ID=putative_sRNA:m", as.integer(x[2]), "_", as.integer(x[3]), ";", sep = ""))

pos_utrs_gr<-subset(red_utrs_gr, strand=="+")
neg_utrs_gr<-subset(red_utrs_gr, strand=="-")

names(pos_utrs_gr) <- apply(as.data.frame(pos_utrs_gr),1, function(x) paste("ID=putative_UTR:p", as.integer(x[2]), "_", as.integer(x[3]), ";", sep = ""))
names(neg_utrs_gr) <- apply(as.data.frame(neg_utrs_gr),1, function(x) paste("ID=putative_UTR:m", as.integer(x[2]), "_", as.integer(x[3]), ";", sep = ""))

# add to major features with strand feature editor (BH function)

# make 'major features file' for each strand (use major_features func from BH)

pos_features<-major_features("PRJNA480455/gff3/original_annotation.gff3", ".", "+", "ncRNA")
neg_features<-major_features("PRJNA480455/gff3/original_annotation.gff3", ".", "-", "ncRNA")

# strand_feature_editor(target_strand, sRNA_IRanges, UTR_IRanges, major_strand_features) 
# this should join utr and srna features together

pos_strand<-strand_feature_editor(target_strand = "+", pos_srnas_gr, pos_utrs_gr, pos_features)
neg_strand<-strand_feature_editor(target_strand = "-", neg_srnas_gr, neg_utrs_gr, neg_features)

# use last step of ffe to create new file

## Creating the final annotation dataframe by combining both strand dataframe and adding missing information like child features from the original GFF3 file.

gff <- read.delim("PRJNA480455/gff3/original_annotation.gff3", header = FALSE, comment.char = "#")

annotation_dataframe <- rbind(gff, pos_strand, neg_strand)

## Remove all the repeating information.

annotation_dataframe <- unique(annotation_dataframe)

## Order the dataframe by feature start coordinates.

annotation_dataframe <- annotation_dataframe[order(annotation_dataframe[,4]),]

## Restore the original header.

f <- readLines("PRJNA480455/gff3/original_annotation.gff3")

header <- c()
i <- 1
while (grepl("#",f[i])==TRUE) {
  f_line <- f[i]
  header <- c(header,f_line)
  i <- i+1
}

# add a line to indicate the origin of the file (single # commas should be ignored by programs)

header <- c(header, "# produced by baerhunter 07_08_23")

## Create the final GFF3 file. 

output_file<-"PRJNA480455/gff3//combined_gffs_07_08_23.gff3"
write.table(header, output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(annotation_dataframe, output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)

### End of modified J. Stiens' code

# Preparing count files for Differential Expression Analysis. Example content of DE_BAM is fixed_SRR7504774/5/6 & fixed_SRR7504780/1/2. (Glucose 24 hr vs Lactate 24hr)
count_features(bam_dir = "PRJNA480455/DE_BAM/", annotation_dir = "PRJNA480455/gff3/", annotation_file = "combined_gffs_07_08_23.gff3", output_dir = "PRJNA480455/Counts/", 
               output_filename = "glu_lac", chromosome_alias_file = "PRJNA480455/aliasfile.txt", is_paired_end = TRUE, strandedness = "stranded")

# The differential expression step - requires a condition (metadata) text file in the format seen on the GitHub. Produces a results.csv file and a DE.res object in R. Taken from Baerhunter.
de.glu_lac <- differential_expression(
  feature_count_file="PRJNA480455/Counts/glu_lac_Counts.csv", 
  metadata_file="PRJNA480455/DE_BAM/glu_lac.txt", 
  cutoff_value=10, multiple_variables=FALSE, main_condition="condition",   
  output_file_name="PRJNA480455/DE_results/glu_lac_results.csv")

# A function to plot create volcano plots based on the results of differential expression:

#' @param de.res A results of differential expression object (obtained from Baerhunter's differential_expression function).
#' @param feature A string indicating which genomic feature to select for: e.g. "UTR", "sRNA", "gene".
#' @param metadata_type A string indicating whether the comparison is made across "time" (6hr vs 24 hr) or "medium" (glucose vs lactate)
#' @param condition A string that goes in front of the genomic feature in the volcano plot title.
#' @param labels A boolean variable that labels the ten most differentially expressed elements in the dataset, provided that they exceed log(2) in differential expression and are statistically significant (padj < 0.05)

de_plot <- function(de.res, feature = "", metadata_type = "", condition = "", labels = TRUE){
    # Reorders the differential expression results by p-value.
    de.res <- de.res[order(de.res$padj),]
    # Adds a column where all elements with adjusted p-value < 0.05 are labelled seperately from the elements that are not statistically significant.
    de.order <- as.data.frame(mutate(as.data.frame(de.res), sig=ifelse(de.res$padj<0.05, "FDR<0.05", "Not Significant"), row.names(rownames(de.res))))
    # Selects only those elements that match the feature variable.
    de.feature <- de.order[grep(feature, rownames(de.order)),]
    # Plots the volcano plot (changes the title format based on metadata_type)
    if (metadata_type == "time") {
      p = ggplot2::ggplot(de.feature, ggplot2::aes(log2FoldChange, -log10(padj))) +
          ggplot2::geom_point(ggplot2::aes(col=sig)) +
          ggplot2::scale_colour_manual(values = c("red", "black")) +
          ggplot2::ggtitle(paste0("Significant features 6hr vs 24hr ", condition, " ", feature))
    }
    if (metadata_type == "medium") {
      p = ggplot2::ggplot(de.feature, ggplot2::aes(log2FoldChange, -log10(padj))) +
          ggplot2::geom_point(ggplot2::aes(col=sig)) +
          ggplot2::scale_colour_manual(values = c("red", "black")) +
          ggplot2::ggtitle(paste0("Significant features ", condition, " ", feature))
    }
    # Collects the names of the ten most differentially expressed genomic elements and assigns them as labels to the appropriate data points.
    if (labels == TRUE) {
      prominent <- de.feature[which(de.feature$padj < 0.05 & abs(de.feature$log2FoldChange) > log2(2)),]
      # Reorders the dataframe so that the most differentially expressed elements are chosen as labels.
      prom2 <- prominent[order(abs(prominent$log2FoldChange), decreasing = TRUE),]
      p + ggrepel::geom_text_repel(data=prom2[1:10,], ggplot2::aes(label=rownames(prom2)[1:10]), force = 100)
    }
    # Occasionally useful to change the function so that the most differentially expressed elements are shown in a dataframe instead of on a graph.
    #return(prom2)
}

# This code is to show the most differentially expressed elements in a dataframe when the above line is not commented out.
glu_lac <- de_plot(de.glu_lac, feature = "UTR", metadata_type = "medium", condition = "glucose vs lactate")
glu_lac[1:10,]
