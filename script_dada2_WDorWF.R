
#Analyse 16s

##### Trimming and clustering sequences #####

#### Dada2 ####

#Installing and reading required packages and libraries

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2")

library(dada2)
library(ShortRead)
library(Biostrings)


path = setwd("C:/Users/lenaf/Documents/R") #Set the working directory path of your project. Your directory should contain your raw sequences in FASTQ format and the cutadapt.exe: https://cutadapt.readthedocs.io/en/stable/installation.html#installation-on-windows

list.files(path)


# Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq and SAMPLENAME_R2.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_RXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_R"), `[`, 1)

#### Identify Primers #### For this step, you should know the sequence of primers you used in your wetlab protocol
FWD = "GTGYCAGCMGCCGCGGTAA"
REV = "GGACTACNVGGGTWTCTAAT"

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

#"pre-filter" the sequences just to remove those with Ns
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = FALSE)

#count the number of times the primers appear
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

#### Remove Primers ####
cutadapt <- "C:/Users/lenaf/Documents/R/cutadapt.exe"  
system2(cutadapt, args = "--version")

#create output filenames for the cutadapt-ed files, and define the parameters we are going to give the cutadapt command.
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 

# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                               "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                               #"--discard-untrimmed", #Discard reads in which no adapter was found
                               fnFs.filtN[i], fnRs.filtN[i])) # input files
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))


#_________________________________________________________#

# Remove reads with more than 'maxLen = x' nucleotides #
fnFs.cutFilt <- file.path(path, "cutFilt", basename(fnFs.cut)) # Put N-filterd files in filtN/ subdirectory
fnRs.cutFilt <- file.path(path, "cutFilt", basename(fnRs.cut))
filterAndTrim(fnFs.cut, fnFs.cutFilt, fnRs.cut, fnRs.cutFilt, maxLen = 390, multithread = FALSE)

path.cutFilt <- file.path(path, "cutFilt")

#_______________________________________________________#

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1.fastq", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format:

sample.names <- sapply(strsplit(basename(fnFs), "_R"), `[`, 1)
sample.names <- sapply(strsplit(basename(sample.names), "01."), `[`, 2)

#### dada2 ####


plotQualityProfile(cutFs[1:2])
plotQualityProfile(cutRs[1:2])

### Based on this plots, you can evaluate where truncate your sequences (truncLen). It should be where 10 nt before the quality starts to decline.

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))


out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, truncLen=c(225,225),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE

head(out)
view(out)

# Learn the Error Rates
# Tool to visualize the frequency of error rate as a function of quality score.
# Necessary for the algorithm - see the paper
errF <- learnErrors(filtFs, multithread=TRUE, MAX_CONSIST = 20)
errR <- learnErrors(filtRs, multithread=TRUE, MAX_CONSIST = 20)



# Dereplication combines all identical sequencing reads into into “unique sequences” with a corresponding 
# “abundance” equal to the number of reads with that unique sequence. Dereplication substantially reduces 
# computation time by eliminating redundant comparisons.

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#We are now ready to apply the core sample inference algorithm to the dereplicated data
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
# dadaFs[[1]]
# dadaRs[[1]]


# Merging paired ends
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
# head(mergers[[1]])

# We can now construct an amplicon sequence variant table (ASV) table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
#table(nchar(getSequences(seqtab)))
#seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(402,427)]
#dim(seqtab2)
#table(nchar(getSequences(seqtab2)))

# Remove chimera
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)


# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.csv(track, file = "Summary.csv")

taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa", minBoot = 80, multithread=TRUE)

## making and writing out standard output files:
# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs_16S <- colnames(seqtab.nochim)
asv_headers_16S <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers_16S[i] <- paste(">ASV_16S", i, sep="_")
}

# fasta:
asv_fasta_16S <- c(rbind(asv_headers_16S, asv_seqs_16S))
write(asv_fasta_16S, "16S_ASVs.fa")

# count table:
asv_tab_16S <- t(seqtab.nochim)
row.names(asv_tab_16S) <- sub(">", "", asv_headers_16S)
write.table(asv_tab_16S, "16S_ASVs_counts.txt", sep="\t", quote=F)

# tax table:
asv_tax_16S <- taxa
row.names(asv_tax_16S) <- sub(">", "", asv_headers_16S)
write.table(asv_tax_16S, "16S_ASVs_taxonomy.txt", sep="\t", quote=F)


#### Clustering OTU 97% ####
library(tibble)
library(dplyr)
library(Biostrings)
library(DECIPHER)

install.packages("installr")
library(installr)
updateR()

install.packages("DECIPHER")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DECIPHER")

nproc <- 12 # set to number of cpus/processors to use for the clustering

asv_sequences <- colnames(seqtab)
sample_names <- rownames(seqtab)
dna <- Biostrings::DNAStringSet(asv_sequences)

installed.packages("DECIPHER")
ls("package:DECIPHER")

## Find clusters of ASVs to form the new OTUs
aln <- DECIPHER::AlignSeqs(dna, processors = nproc)
d <- DECIPHER::DistanceMatrix(aln, processors = nproc)
clusters <- DECIPHER::TreeLine(
  myDistMatrix=d,
  method = "complete",
  cutoff = 0.03, # use `cutoff = 0.03` for a 97% OTU
  type = "clusters",
  processors = nproc)

nrow(clusters)
nrow(asv_sequences)


## Use dplyr to merge the columns of the seqtab matrix for ASVs in the same OTU
# prep by adding sequences to the `clusters` data frame
clusters <- clusters %>%
  add_column(sequence = asv_sequences)
merged_seqtab <- seqtab %>%
  # setup: turn seqtab into a tibble with rows = ASVs and columns = samples
  t %>%
  as_tibble(rownames = "sequence") %>%
  # add the cluster information
  left_join(clusters, by = "sequence") %>%
  # merge ASVs in the same cluster, summing abundances within samples
  group_by(cluster) %>%
  summarize_at(vars(-sequence), sum) %>%
  # Set new taxa names to OTU<cluster #> 
  mutate(cluster = paste0("OTU", cluster)) %>%
  # convert back to a matrix in the original orientation
  column_to_rownames("cluster") %>%
  as("matrix") %>%
  t

  # count table:

write.table(merged_seqtab, "Mito_OTUs_counts.txt", sep="\t", quote=F)

#### Match a reference sequence from ASVs table to the corresponding OTU 95% ####
#library(dplyr)
newdf<-as.data.frame(cbind.data.frame(rownames(clusters),clusters$cluster,clusters$sequence))
colnames(newdf)<-c("ASV","C.ASV","sequence")

newdf.sort<-newdf %>% arrange(C.ASV) %>% distinct(C.ASV,.keep_all = TRUE)

  #FASTA file
fas<-paste(">OTU",seq(1:nrow(newdf.sort)),"\n",newdf.sort$sequence,sep="")
write.table(fas,file="OTU.fas", row.names=FALSE, col.names=FALSE,quote = FALSE)


#### OTU Taxonomic assignation ####

merged_seqtab2 <- merged_seqtab
colnames(merged_seqtab2) <-newdf.sort$sequence #Change Column name with Sequences of corresponding OTU

taxaOTU <- assignTaxonomy(merged_seqtab2, "C:/Users/lenaf/Documents/R/silva_nr99_v138.1_train_set.fa", multithread = TRUE)

  # tax table:
OTU_headers <- vector(dim(merged_seqtab2)[2], mode="character")
for (i in 1:dim(merged_seqtab2)[2]) {
  OTU_headers[i] <- paste("Mito OTU", i, sep="")
}
OTU_tax <- taxaOTU
row.names(OTU_tax) <- OTU_headers
write.table(OTU_tax, "Mito_OTUs_taxonomy.txt", sep="\t", quote=F)

#At the end of this analyse, you should have obtained:

#Mito_OTUs_taxonomy.txt: This file contains the taxonomic assignments for the mitochondrial OTUs.
#OTU.fas: This file contains the representative sequences for each OTU.
#Mito_OTUs_counts.txt: This file contains the OTU abundance table for the mitochondrial OTUs.
#16S_ASVs_taxonomy.txt: This file contains the taxonomic assignments for the 16S rRNA gene ASVs.
#16S_ASVs_counts.txt: This file contains the ASV abundance table for the 16S rRNA gene.
#16S_ASVs.fa: This file contains the representative sequences for each 16S rRNA gene ASV.
#Summary.csv: This file contains summary statistics for the DADA2 analysis.
#The folders cutadapt, cultFilt, and filtN that contain intermediate files generated during the quality control and filtering steps of the DADA2 analysis.



############################################## Diversity analysis with Phyloseq ####

R.Version()$version.string #The version 4.2.2 of R is not compatible with phyloseq. Use the version 4.0.1 instead (windows 64bit)

path = setwd("C:/Users/lenaf/Documents/R") #Accesing the working directory

update.packages(checkBuilt = TRUE, ask = FALSE)

# Load ASV abundance table and taxonomic assignment files
asv_table <- read.table("16S_ASVs_counts.txt", header = TRUE, row.names = 1, sep = "\t")
asv_taxonomy <- read.table("16S_ASVs_taxonomy.txt", header = TRUE, row.names = 1, sep = "\t")

#Installing and reading required packages and libraries

install.packages("BiocManager")
BiocManager::install("phyloseq")
update.packages()

library(BiocManager)
library(phyloseq)
install.packages('ggplot2')
library(ggplot2)
install.packages('vegan')
library(vegan)
install.packages('dplyr')
library(dplyr)
library(scales)
library(grid)
install.packages('reshape2')
library(reshape2)
install.packages('cowplot')
library(cowplot)
install.packages('BiocManager')

options(install.packages.check.source = "no")
BiocManager::install("phyloseq")
library(rvest)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")
install.packages("Bioconductor")
library(BiocManager)
source("https://raw.githubusercontent.com/joey711/phyloseq/master/inst/scripts/installer.R",
       local = TRUE)

library(ggplot2)
theme_set(theme_bw())


#You will use ASVs files to create your OTU_table
# Load ASV abundance table and taxonomic assignment files
# Convert taxonomy table to matrix and create phyloseq object
tax_mat <- as.matrix(asv_taxonomy)
ps <- phyloseq(otu_table(asv_table, taxa_are_rows = TRUE), tax_table(tax_mat))

# Extract the sample names from the existing names (leaving just the sampleID)
asv_table_sample_names <- sub("^[^.]*\\.", "", colnames(asv_table))

# Update the column names in the ASV table
colnames(asv_table) <- asv_table_sample_names

# Check if taxa names match
all(taxa_names(asv_table) == taxa_names(asv_taxonomy)) # If this returns FALSE, there is a mismatch in taxa names.

# Check order of taxa
all.equal(taxa_names(asv_table), rownames(otu_table)) # If this returns FALSE, there is a mismatch in the order of taxa.

# Convert taxonomy table to matrix
tax_mat <- as.matrix(asv_taxonomy)

# Create phyloseq object
ps <- phyloseq(otu_table(asv_table, taxa_are_rows = TRUE), tax_table(tax_mat))

# Collapse ASVs at the same taxonomic level
ps_glom <- tax_glom(ps, "Class")

# Get the OTU table
otu_table <- as.matrix(otu_table(ps_glom))

# Save the OTU table
write.table(otu_table, file = "ASV_OTU_table.txt", sep = "\t")

# Load the metadata file
Metadata <- read.delim("Summary2.txt", header = TRUE, row.names = 1) #Summary2.txt is a file that I manually created based on Summary.csv. I added metadata information (data of sample collection, collection site etc) manually and than saved the file in .txt. You should create a metadata file for your samples.

# Extract the sample names from the existing names (leaving just the sampleID)
Metadata_sample_names <- sub("^[^.]*\\.", "", rownames(Metadata))

# Update the row names in the metadata file
rownames(Metadata) <- Metadata_sample_names

# Check if the sample names match
identical(asv_table_sample_names, Metadata_sample_names)

# Use match() to find the indices of the metadata samples in the ASV table
sample_indices <- match(Metadata_sample_names, asv_table_sample_names)

#You can associate your phyloseq object with your metadata file. 
#First, you need to make sure that the sample names in your metadata
#file match the sample names in your phyloseq object. You can check this 
#by running:

all(sample_names(Metadata) %in% sample_names(ps)) #If the result is TRUE, then the sample names match

setdiff(sample_names(Metadata), sample_names(ps))
setdiff(sample_names(ps), sample_names(Metadata))

all(sample_names(Metadata) %in% sample_names(ps))
setdiff(sample_names(Metadata), sample_names(ps))
setdiff(sample_names(ps), sample_names(Metadata))

#Once you have confirmed that the sample names match, you can add the metadata to your phyloseq object like this:

ps <- merge_phyloseq(ps, Metadata)


# Get the sample names from the phyloseq object
sample_names_ps <- sample_names(ps)

# Get the sample names from the metadata file
sample_names_metadata <- rownames(Metadata)

# Compare the sample names
identical(sample_names_ps, sample_names_metadata)

#########Yes, if the output of identical(sample_names(ps), 
#rownames(metadata)) is TRUE, it means that the sample names in
#your phyloseq object match the sample names in your metadata, 
#and you have successfully created the phyloseq object with metadata.
#Once you have created your phyloseq object with metadata, you can 
#perform a wide range of analysis on your microbiome data. 

#Verify the phyloseq version
packageVersion("phyloseq")

#verify the functions on Phyloseq
ls("package:phyloseq")

################### Calculate the total number of reads per sample
sample_sums <- as.data.frame(sample_sums(ps))

# Rename the column to "read_count"
colnames(sample_sums) <- "read_count"

# Create a histogram of the read counts
ggplot(sample_sums, aes(x = read_count)) +
  geom_histogram(binwidth = 10000, fill = "#87CEFA", color = "black") +
  labs(title = "Distribution of read counts per sample",
       x = "Read count",
       y = "Frequency")

# calculate the minimum, average and maximum sample read counts
summary(sample_sums$read_count)

################Calculate alpha diversity measures (all samples):
alpha_div <- estimate_richness(ps)

# Save alpha diversity results as a CSV file
write.csv(alpha_div, file = "alpha_diversity.csv", row.names = TRUE)

#################### Ordination analysis - by Bray Curtis distance #####################

ordination <- ordinate(ps, method = "PCoA", distance = "bray")

plot_ordination(ps, ordination, color = "SampleType")

# Install the lubridate package
install.packages("lubridate")

# Load the lubridate package
library(lubridate)
library(ggplot2)
library(phyloseq)


# Add a new column for season
Metadata$season <- ifelse(month(Metadata$Date) %in% c(12, 1, 2), "Winter",
                          ifelse(month(Metadata$Date) %in% c(3, 4, 5), "Spring",
                                 ifelse(month(Metadata$Date) %in% c(6, 7, 8), "Summer", "Fall")))


# Assign the sample metadata to the sam_data slot of the ps object
sample_data(ps) <- Metadata

#######Ordination by wastewater station

# Assign the factor variable to the SampleType column in the sample metadata
sample_data(ps)$SampleType <- factor(sample_data(ps)$POINT)

# Ordination analysis
ordination <- ordinate(ps, method = "PCoA", distance = "bray")

# Plot ordination with points colored by SampleType (POINT variable)
plot_ordination(ps, ordination, color = "SampleType") +
  geom_point(size = 3) +
  geom_text(aes(label = sample_names(ps)), size = 3, nudge_y = 0.02) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic()

###### Ordination grouped by season

# Assign the factor variable to the SampleType column in the sample metadata
sample_data(ps)$SampleType <- factor(sample_data(ps)$season)

# Ordination analysis
ordination <- ordinate(ps, method = "PCoA", distance = "bray")

# Plot ordination with points colored by SampleType (season variable)
plot_ordination(ps, ordination, color = "SampleType") +
  geom_point(size = 3) +
  geom_text(aes(label = sample_names(ps)), size = 3, nudge_y = 0.02) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic()

###### Ordination grouped by concentration method

# Assign the factor variable to the SampleType column in the sample metadata
sample_data(ps)$SampleType <- factor(sample_data(ps)$MET)

# Ordination analysis
ordination <- ordinate(ps, method = "PCoA", distance = "bray")

# Plot ordination with points colored by SampleType (MET variable)
plot_ordination(ps, ordination, color = "SampleType") +
  geom_point(size = 3) +
  geom_text(aes(label = sample_names(ps)), size = 3, nudge_y = 0.02) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic()

###### Ordination grouped by Extraction kit

# Assign the factor variable to the SampleType column in the sample metadata
sample_data(ps)$SampleType <- factor(sample_data(ps)$KIT)

# Ordination analysis
ordination <- ordinate(ps, method = "PCoA", distance = "bray")

# Plot ordination with points colored by SampleType (KIT variable)
plot_ordination(ps, ordination, color = "SampleType") +
  geom_point(size = 3) +
  geom_text(aes(label = sample_names(ps)), size = 3, nudge_y = 0.02) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic()

###### Ordination grouped by centrifugation speed

# Assign the factor variable to the SampleType column in the sample metadata
sample_data(ps)$SampleType <- factor(sample_data(ps)$VIT)

# Ordination analysis
ordination <- ordinate(ps, method = "PCoA", distance = "bray")

# Plot ordination with points colored by SampleType (VIT variable)
plot_ordination(ps, ordination, color = "SampleType") +
  geom_point(size = 3) +
  geom_text(aes(label = sample_names(ps)), size = 3, nudge_y = 0.02) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic()

###### Ordination grouped by Volume

# Assign the factor variable to the SampleType column in the sample metadata
sample_data(ps)$SampleType <- factor(sample_data(ps)$VOL)

# Ordination analysis
ordination <- ordinate(ps, method = "PCoA", distance = "bray")

# Plot ordination with points colored by SampleType (POINT variable)
plot_ordination(ps, ordination, color = "SampleType") +
  geom_point(size = 3) +
  geom_text(aes(label = sample_names(ps)), size = 3, nudge_y = 0.02) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic()

######################################## Essai volume

# Subset the phyloseq object for samples collected on December 7, 2021
ps_subset <- subset_samples(ps, sample_data(ps)$Date == "2021-12-07")

# Assign the factor variable to the POINT column in the sample metadata
sample_data(ps_subset)$SampleType <- factor(sample_data(ps_subset)$POINT)

# Perform ordination analysis
ordination <- ordinate(ps_subset, method = "PCoA", distance = "bray")

# Plot the ordination with points colored by SampleType (POINT variable)
plot_ordination(ps_subset, ordination, color = "SampleType") +
  geom_point(size = 3) +
  geom_text(aes(label = sample_names(ps_subset)), size = 3, vjust = 1.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic()

################ Comparisons between samples pairs 
# Ensure we are using the subset of the data from December 7, 2021
ps_subset <- subset_samples(ps, sample_data(ps)$Date == "2021-12-07")

sample_names(ps_subset)


# Calculate the Bray distance matrix
bray_dist_matrix <- phyloseq::distance(ps_subset, method = "bray")

# Convert the distance matrix into a full matrix
bray_dist_full <- as.matrix(bray_dist_matrix)

# Extract the dissimilarities for the specific pairs of samples using the sample names from ps_subset
dissimilarity_34_35 <- bray_dist_full["34", "35"] #PELLET 50 VS PELLET 125
dissimilarity_34_38 <- bray_dist_full["34", "38"] #PELLET 50 VS FILTER 50
dissimilarity_35_37 <- bray_dist_full["35", "37"] #PELLET 125 VS FILTER 250
dissimilarity_38_37 <- bray_dist_full["38", "37"] #FILTER 50 VS FILTER 250

dissimilarity_39_41 <- bray_dist_full["39", "41"] #PELLET 50 VS PELLET 125
dissimilarity_39_43 <- bray_dist_full["39", "43"] #PELLET 50 VS FILTER 50
dissimilarity_41_42 <- bray_dist_full["41", "42"] #PELLET 125 VS FILTER 250
dissimilarity_43_42 <- bray_dist_full["43", "42"] #FILTER 50 VS FILTER 250

dissimilarity_44_46 <- bray_dist_full["44", "46"] #PELLET 50 VS FILTER 125
dissimilarity_44_47 <- bray_dist_full["44", "47"] #PELLET 50 VS FILTER 50
dissimilarity_46_47 <- bray_dist_full["46", "47"] #FILTER 125 VS FILTER 50

# Create a dataframe to nicely format the extracted dissimilarities
dissimilarity_values <- data.frame(
  Sample_Pair = c("34-35", "34-38", "35-37", "38-37", 
                  "39-41", "39-43", "41-42", "43-42", 
                  "44-46", "44-47", "46-47"),
  Bray_Dissimilarity = c(dissimilarity_34_35, dissimilarity_34_38, 
                            dissimilarity_35_37, dissimilarity_38_37, 
                            dissimilarity_39_41, dissimilarity_39_43, 
                            dissimilarity_41_42, dissimilarity_43_42, 
                            dissimilarity_44_46, dissimilarity_44_47, 
                            dissimilarity_46_47)
)
# Print the dissimilarity values
print(dissimilarity_values)

################################## Essai Kit

# Subset the phyloseq object for samples collected on 2022-01-20
ps_subset <- subset_samples(ps, sample_data(ps)$Date == "2022-01-20")

# Assign the factor variable to the POINT column in the sample metadata
sample_data(ps_subset)$SampleType <- factor(sample_data(ps_subset)$POINT)

# Perform ordination analysis
ordination <- ordinate(ps_subset, method = "PCoA", distance = "bray")

# Plot the ordination with points colored by SampleType (POINT variable)
plot_ordination(ps_subset, ordination, color = "SampleType") +
  geom_point(size = 3) +
  geom_text(aes(label = sample_names(ps_subset)), size = 3, vjust = 1.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic()

################ Comparisons between samples pairs 
# Ensure we are using the subset of the data from 2022-01-20
ps_subset <- subset_samples(ps, sample_data(ps)$Date == "2022-01-20")

sample_names(ps_subset)


# Calculate the Bray distance matrix
bray_dist_matrix <- phyloseq::distance(ps_subset, method = "bray")

# Convert the distance matrix into a full matrix
bray_dist_full <- as.matrix(bray_dist_matrix)

# Extract the dissimilarities for the specific pairs of samples using the sample names from ps_subset
dissimilarity_48_49 <- bray_dist_full["48", "49"] #PELLET PowerLyzer VS PELLET PowerSoil Pro
dissimilarity_48_51 <- bray_dist_full["48", "51"] #PELLET Powerlyzer VS FILTER Powerlyzer
dissimilarity_50_51 <- bray_dist_full["50", "51"] #FILTER PowersoilPro VS FILTER PowerLyzer
dissimilarity_50_49 <- bray_dist_full["50", "49"] #FILTER PowerSoil Pro VS PELLET PowerSoil Pro

dissimilarity_52_53 <- bray_dist_full["52", "53"] #PELLET PowerLyzer VS PELLET PowerSoil Pro
dissimilarity_52_55 <- bray_dist_full["52", "55"] #PELLET Powerlyzer VS FILTER Powerlyzer
dissimilarity_54_55 <- bray_dist_full["54", "55"] #FILTER PowersoilPro VS FILTER PowerLyzer
dissimilarity_54_53 <- bray_dist_full["54", "53"] #FILTER PowerSoil Pro VS PELLET PowerSoil Pro

dissimilarity_56_57 <- bray_dist_full["56", "57"] #PELLET PowerLyzer VS PELLET PowerSoil Pro
dissimilarity_58_57 <- bray_dist_full["58", "57"] #FILTER PowerSoil Pro VS PELLET PowerSoil Pro
dissimilarity_58_56 <- bray_dist_full["58", "56"] #FILTER PowerSoil Pro VS PELLET PowerLyzer 


# Create a dataframe to nicely format the extracted dissimilarities
dissimilarity_values <- data.frame(
  Sample_Pair = c("48-49", "48-51", "50-51", "50-49", 
                  "52-53", "52-55", "54-55", "54-53", 
                  "56-57", "58-57", "58-56"),
  Bray_Dissimilarity = c(dissimilarity_48_49, dissimilarity_48_51, 
                         dissimilarity_50_51, dissimilarity_50_49, 
                         dissimilarity_52_53, dissimilarity_52_55, 
                         dissimilarity_54_55, dissimilarity_54_53, 
                         dissimilarity_56_57, dissimilarity_58_57, 
                         dissimilarity_58_56)
)
# Print the dissimilarity values
print(dissimilarity_values)

############################  Essai Vitesse

# Subset the phyloseq object for samples collected on 2022-02-17
ps_subset <- subset_samples(ps, sample_data(ps)$Date == "2022-02-17")

# Assign the factor variable to the POINT column in the sample metadata
sample_data(ps_subset)$SampleType <- factor(sample_data(ps_subset)$POINT)

# Perform ordination analysis
ordination <- ordinate(ps_subset, method = "PCoA", distance = "bray")

# Plot the ordination with points colored by SampleType (POINT variable)
plot_ordination(ps_subset, ordination, color = "SampleType") +
  geom_point(size = 3) +
  geom_text(aes(label = sample_names(ps_subset)), size = 3, vjust = 1.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic()

################ Comparisons between samples pairs 
# Ensure we are using the subset of the data from 2022-02-17
ps_subset <- subset_samples(ps, sample_data(ps)$Date == "2022-02-17")

sample_names(ps_subset)


# Calculate the Bray distance matrix
bray_dist_matrix <- phyloseq::distance(ps_subset, method = "bray")

# Convert the distance matrix into a full matrix
bray_dist_full <- as.matrix(bray_dist_matrix)

# Extract the dissimilarities for the specific pairs of samples using the sample names from ps_subset
dissimilarity_60_62 <- bray_dist_full["60", "62"] #PELLET 3400 VS PELLET 15 000

dissimilarity_64_66 <- bray_dist_full["64", "66"] #PELLET 3400 VS PELLET 15 000

dissimilarity_68_70 <- bray_dist_full["68", "70"] #PELLET 3400 VS PELLET 15 000

# Create a dataframe to nicely format the extracted dissimilarities
dissimilarity_values <- data.frame(
  Sample_Pair = c("60-62", "64-66", "68-70"),
  Bray_Dissimilarity = c(dissimilarity_60_62, dissimilarity_64_66, dissimilarity_68_70)
)

# Print the dissimilarity values
print(dissimilarity_values)


######################## Essai Filtre, Culot, Culot+Filtre

# Subset the phyloseq object for samples collected on 2022-04-21
ps_subset <- subset_samples(ps, sample_data(ps)$Date == "2022-04-21")

# Assign the factor variable to the POINT column in the sample metadata
sample_data(ps_subset)$SampleType <- factor(sample_data(ps_subset)$POINT)

# Perform ordination analysis
ordination <- ordinate(ps_subset, method = "PCoA", distance = "bray")

# Plot the ordination with points colored by SampleType (POINT variable)
plot_ordination(ps_subset, ordination, color = "SampleType") +
  geom_point(size = 3) +
  geom_text(aes(label = sample_names(ps_subset)), size = 3, vjust = 1.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic()

################ Comparisons between samples pairs 
# Ensure we are using the subset of the data from 2022-04-21
ps_subset <- subset_samples(ps, sample_data(ps)$Date == "2022-04-21")

sample_names(ps_subset)

# Calculate the Bray distance matrix
bray_dist_matrix <- phyloseq::distance(ps_subset, method = "bray")

# Convert the distance matrix into a full matrix
bray_dist_full <- as.matrix(bray_dist_matrix)

# Extract the dissimilarities for the specific pairs of samples using the sample names from ps_subset
dissimilarity_87_89 <- bray_dist_full["87", "89"] #FILTER VS PELLET
dissimilarity_87_88 <- bray_dist_full["87", "88"] #FILTER VS PELLET+FILTER
dissimilarity_89_88 <- bray_dist_full["89", "88"] #PELLET VS PELLET+FILTER

dissimilarity_90_92 <- bray_dist_full["90", "92"] #FILTER VS PELLET
dissimilarity_90_91 <- bray_dist_full["90", "91"] #FILTER VS PELLET+FILTER
dissimilarity_92_91 <- bray_dist_full["92", "91"] #PELLET VS PELLET+FILTER
dissimilarity_92_93 <- bray_dist_full["92", "93"] #PELLET VS PELLET+MGCL2

dissimilarity_94_96 <- bray_dist_full["94", "96"] #FILTER VS PELLET
dissimilarity_94_95 <- bray_dist_full["94", "95"] #FILTER VS PELLET+FILTER
dissimilarity_96_95 <- bray_dist_full["96", "95"] #PELLET VS PELLET+FILTER

# Create a dataframe to nicely format the extracted dissimilarities
dissimilarity_values <- data.frame(
  Sample_Pair = c("87-89", "87-88", "89-88", "90-92", "90-91", "92-91", "92-93", "94-96", "94-95", "96-95"),
  Bray_Dissimilarity = c(dissimilarity_87_89, dissimilarity_87_88, dissimilarity_89_88, 
                         dissimilarity_90_92, dissimilarity_90_91, dissimilarity_92_91, 
                         dissimilarity_92_93, dissimilarity_94_96, dissimilarity_94_95, dissimilarity_96_95)
)

# Print the dissimilarity values
print(dissimilarity_values)

########################## Ordination with Jaccard ##############

ordination <- ordinate(ps, method = "PCoA", distance = "jaccard")

plot_ordination(ps, ordination, color = "SampleType")

# Install the lubridate package
install.packages("lubridate")

# Load the lubridate package
library(lubridate)
library(ggplot2)
library(phyloseq)


# Add a new column for season
Metadata$season <- ifelse(month(Metadata$Date) %in% c(12, 1, 2), "Winter",
                          ifelse(month(Metadata$Date) %in% c(3, 4, 5), "Spring",
                                 ifelse(month(Metadata$Date) %in% c(6, 7, 8), "Summer", "Fall")))


# Assign the sample metadata to the sam_data slot of the ps object
sample_data(ps) <- Metadata

#######Ordination by wastewater station

# Assign the factor variable to the SampleType column in the sample metadata
sample_data(ps)$SampleType <- factor(sample_data(ps)$POINT)

# Ordination analysis
ordination <- ordinate(ps, method = "PCoA", distance = "jaccard")

# Plot ordination with points colored by SampleType (POINT variable)
plot_ordination(ps, ordination, color = "SampleType") +
  geom_point(size = 3) +
  geom_text(aes(label = sample_names(ps)), size = 3, nudge_y = 0.02) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic()

###### Ordination grouped by season

# Assign the factor variable to the SampleType column in the sample metadata
sample_data(ps)$SampleType <- factor(sample_data(ps)$season)

# Ordination analysis
ordination <- ordinate(ps, method = "PCoA", distance = "jaccard")

# Plot ordination with points colored by SampleType (season variable)
plot_ordination(ps, ordination, color = "SampleType") +
  geom_point(size = 3) +
  geom_text(aes(label = sample_names(ps)), size = 3, nudge_y = 0.02) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic()

###### Ordination grouped by concentration method

# Assign the factor variable to the SampleType column in the sample metadata
sample_data(ps)$SampleType <- factor(sample_data(ps)$MET)

# Ordination analysis
ordination <- ordinate(ps, method = "PCoA", distance = "jaccard")

# Plot ordination with points colored by SampleType (MET variable)
plot_ordination(ps, ordination, color = "SampleType") +
  geom_point(size = 3) +
  geom_text(aes(label = sample_names(ps)), size = 3, nudge_y = 0.02) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic()

###### Ordination grouped by Extraction kit

# Assign the factor variable to the SampleType column in the sample metadata
sample_data(ps)$SampleType <- factor(sample_data(ps)$KIT)

# Ordination analysis
ordination <- ordinate(ps, method = "PCoA", distance = "jaccard")

# Plot ordination with points colored by SampleType (KIT variable)
plot_ordination(ps, ordination, color = "SampleType") +
  geom_point(size = 3) +
  geom_text(aes(label = sample_names(ps)), size = 3, nudge_y = 0.02) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic()

###### Ordination grouped by centrifugation speed

# Assign the factor variable to the SampleType column in the sample metadata
sample_data(ps)$SampleType <- factor(sample_data(ps)$VIT)

# Ordination analysis
ordination <- ordinate(ps, method = "PCoA", distance = "jaccard")

# Plot ordination with points colored by SampleType (VIT variable)
plot_ordination(ps, ordination, color = "SampleType") +
  geom_point(size = 3) +
  geom_text(aes(label = sample_names(ps)), size = 3, nudge_y = 0.02) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic()

###### Ordination grouped by Volume

# Assign the factor variable to the SampleType column in the sample metadata
sample_data(ps)$SampleType <- factor(sample_data(ps)$VOL)

# Ordination analysis
ordination <- ordinate(ps, method = "PCoA", distance = "jaccard")

# Plot ordination with points colored by SampleType (POINT variable)
plot_ordination(ps, ordination, color = "SampleType") +
  geom_point(size = 3) +
  geom_text(aes(label = sample_names(ps)), size = 3, nudge_y = 0.02) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic()

######################################## Essai volume

# Subset the phyloseq object for samples collected on December 7, 2021
ps_subset <- subset_samples(ps, sample_data(ps)$Date == "2021-12-07")

# Assign the factor variable to the POINT column in the sample metadata
sample_data(ps_subset)$SampleType <- factor(sample_data(ps_subset)$POINT)

# Perform ordination analysis
ordination <- ordinate(ps_subset, method = "PCoA", distance = "jaccard")

# Plot the ordination with points colored by SampleType (POINT variable)
plot_ordination(ps_subset, ordination, color = "SampleType") +
  geom_point(size = 3) +
  geom_text(aes(label = sample_names(ps_subset)), size = 3, vjust = 1.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic()

library(phyloseq)
library(vegan)


################ Comparisons between samples pairs 
# Ensure we are using the subset of the data from December 7, 2021
ps_subset <- subset_samples(ps, sample_data(ps)$Date == "2021-12-07")

sample_names(ps_subset)


# Calculate the Jaccard distance matrix
jaccard_dist_matrix <- phyloseq::distance(ps_subset, method = "jaccard")

# Convert the distance matrix into a full matrix
jaccard_dist_full <- as.matrix(jaccard_dist_matrix)

# Extract the dissimilarities for the specific pairs of samples using the sample names from ps_subset
dissimilarity_34_35 <- jaccard_dist_full["34", "35"] #PELLET 50 VS PELLET 125
dissimilarity_34_38 <- jaccard_dist_full["34", "38"] #PELLET 50 VS FILTER 50
dissimilarity_35_37 <- jaccard_dist_full["35", "37"] #PELLET 125 VS FILTER 250
dissimilarity_38_37 <- jaccard_dist_full["38", "37"] #FILTER 50 VS FILTER 250

dissimilarity_39_41 <- jaccard_dist_full["39", "41"] #PELLET 50 VS PELLET 125
dissimilarity_39_43 <- jaccard_dist_full["39", "43"] #PELLET 50 VS FILTER 50
dissimilarity_41_42 <- jaccard_dist_full["41", "42"] #PELLET 125 VS FILTER 250
dissimilarity_43_42 <- jaccard_dist_full["43", "42"] #FILTER 50 VS FILTER 250

dissimilarity_44_46 <- jaccard_dist_full["44", "46"] #PELLET 50 VS FILTER 125
dissimilarity_44_47 <- jaccard_dist_full["44", "47"] #PELLET 50 VS FILTER 50
dissimilarity_46_47 <- jaccard_dist_full["46", "47"] #FILTER 125 VS FILTER 50

# Create a dataframe to nicely format the extracted dissimilarities
dissimilarity_values <- data.frame(
  Sample_Pair = c("34-35", "34-38", "35-37", "38-37", 
                  "39-41", "39-43", "41-42", "43-42", 
                  "44-46", "44-47", "46-47"),
  Jaccard_Dissimilarity = c(dissimilarity_34_35, dissimilarity_34_38, 
                            dissimilarity_35_37, dissimilarity_38_37, 
                            dissimilarity_39_41, dissimilarity_39_43, 
                            dissimilarity_41_42, dissimilarity_43_42, 
                            dissimilarity_44_46, dissimilarity_44_47, 
                            dissimilarity_46_47)
)
# Print the dissimilarity values
print(dissimilarity_values)

################################## Essai Kit

# Subset the phyloseq object for samples collected on December 7, 2021
ps_subset <- subset_samples(ps, sample_data(ps)$Date == "2022-01-20")

# Assign the factor variable to the POINT column in the sample metadata
sample_data(ps_subset)$SampleType <- factor(sample_data(ps_subset)$POINT)

# Perform ordination analysis
ordination <- ordinate(ps_subset, method = "PCoA", distance = "jaccard")

# Plot the ordination with points colored by SampleType (POINT variable)
plot_ordination(ps_subset, ordination, color = "SampleType") +
  geom_point(size = 3) +
  geom_text(aes(label = sample_names(ps_subset)), size = 3, vjust = 1.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic()

################ Comparisons between samples pairs 
# Ensure we are using the subset of the data from December 7, 2021
ps_subset <- subset_samples(ps, sample_data(ps)$Date == "2022-01-20")

sample_names(ps_subset)


# Calculate the jaccard distance matrix
jaccard_dist_matrix <- phyloseq::distance(ps_subset, method = "jaccard")

# Convert the distance matrix into a full matrix
jaccard_dist_full <- as.matrix(jaccard_dist_matrix)

# Extract the dissimilarities for the specific pairs of samples using the sample names from ps_subset
dissimilarity_48_49 <- jaccard_dist_full["48", "49"] #PELLET PowerLyzer VS PELLET PowerSoil Pro
dissimilarity_48_51 <- jaccard_dist_full["48", "51"] #PELLET Powerlyzer VS FILTER Powerlyzer
dissimilarity_50_51 <- jaccard_dist_full["50", "51"] #FILTER PowersoilPro VS FILTER PowerLyzer
dissimilarity_50_49 <- jaccard_dist_full["50", "49"] #FILTER PowerSoil Pro VS PELLET PowerSoil Pro

dissimilarity_52_53 <- jaccard_dist_full["52", "53"] #PELLET PowerLyzer VS PELLET PowerSoil Pro
dissimilarity_52_55 <- jaccard_dist_full["52", "55"] #PELLET Powerlyzer VS FILTER Powerlyzer
dissimilarity_54_55 <- jaccard_dist_full["54", "55"] #FILTER PowersoilPro VS FILTER PowerLyzer
dissimilarity_54_53 <- jaccard_dist_full["54", "53"] #FILTER PowerSoil Pro VS PELLET PowerSoil Pro

dissimilarity_56_57 <- jaccard_dist_full["56", "57"] #PELLET PowerLyzer VS PELLET PowerSoil Pro
dissimilarity_58_57 <- jaccard_dist_full["58", "57"] #FILTER PowerSoil Pro VS PELLET PowerSoil Pro
dissimilarity_58_56 <- jaccard_dist_full["58", "56"] #FILTER PowerSoil Pro VS PELLET PowerLyzer 


# Create a dataframe to nicely format the extracted dissimilarities
dissimilarity_values <- data.frame(
  Sample_Pair = c("48-49", "48-51", "50-51", "50-49", 
                  "52-53", "52-55", "54-55", "54-53", 
                  "56-57", "58-57", "58-56"),
  Jaccard_Dissimilarity = c(dissimilarity_48_49, dissimilarity_48_51, 
                         dissimilarity_50_51, dissimilarity_50_49, 
                         dissimilarity_52_53, dissimilarity_52_55, 
                         dissimilarity_54_55, dissimilarity_54_53, 
                         dissimilarity_56_57, dissimilarity_58_57, 
                         dissimilarity_58_56)
)
# Print the dissimilarity values
print(dissimilarity_values)


############################  Essai Vitesse

# Subset the phyloseq object for samples collected on December 7, 2021
ps_subset <- subset_samples(ps, sample_data(ps)$Date == "2022-02-17")

# Assign the factor variable to the POINT column in the sample metadata
sample_data(ps_subset)$SampleType <- factor(sample_data(ps_subset)$POINT)

# Perform ordination analysis
ordination <- ordinate(ps_subset, method = "PCoA", distance = "jaccard")

# Plot the ordination with points colored by SampleType (POINT variable)
plot_ordination(ps_subset, ordination, color = "SampleType") +
  geom_point(size = 3) +
  geom_text(aes(label = sample_names(ps_subset)), size = 3, vjust = 1.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic()


# Ensure we are using the subset of the data from 2022-02-17
ps_subset <- subset_samples(ps, sample_data(ps)$Date == "2022-02-17")

sample_names(ps_subset)


# Calculate the jaccard distance matrix
jaccard_dist_matrix <- phyloseq::distance(ps_subset, method = "jaccard")

# Convert the distance matrix into a full matrix
jaccard_dist_full <- as.matrix(jaccard_dist_matrix)

# Extract the dissimilarities for the specific pairs of samples using the sample names from ps_subset
dissimilarity_60_62 <- jaccard_dist_full["60", "62"] #PELLET 3400 VS PELLET 15 000

dissimilarity_64_66 <- jaccard_dist_full["64", "66"] #PELLET 3400 VS PELLET 15 000

dissimilarity_68_70 <- jaccard_dist_full["68", "70"] #PELLET 3400 VS PELLET 15 000

# Create a dataframe to nicely format the extracted dissimilarities
dissimilarity_values <- data.frame(
  Sample_Pair = c("60-62", "64-66", "68-70"),
  Jaccard_Dissimilarity = c(dissimilarity_60_62, dissimilarity_64_66, dissimilarity_68_70)
)

# Print the dissimilarity values
print(dissimilarity_values)

######################## Essai Filtre, Culot, Culot+Filtre

# Subset the phyloseq object for samples collected on December 7, 2021
ps_subset <- subset_samples(ps, sample_data(ps)$Date == "2022-04-21")

# Assign the factor variable to the POINT column in the sample metadata
sample_data(ps_subset)$SampleType <- factor(sample_data(ps_subset)$POINT)

# Perform ordination analysis
ordination <- ordinate(ps_subset, method = "PCoA", distance = "jaccard")

# Plot the ordination with points colored by SampleType (POINT variable)
plot_ordination(ps_subset, ordination, color = "SampleType") +
  geom_point(size = 3) +
  geom_text(aes(label = sample_names(ps_subset)), size = 3, vjust = 1.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic()

################ Comparisons between samples pairs 
# Ensure we are using the subset of the data from 2022-04-21
ps_subset <- subset_samples(ps, sample_data(ps)$Date == "2022-04-21")

sample_names(ps_subset)

# Calculate the jaccard distance matrix
jaccard_dist_matrix <- phyloseq::distance(ps_subset, method = "jaccard")

# Convert the distance matrix into a full matrix
jaccard_dist_full <- as.matrix(jaccard_dist_matrix)

# Extract the dissimilarities for the specific pairs of samples using the sample names from ps_subset
dissimilarity_87_89 <- jaccard_dist_full["87", "89"] #FILTER VS PELLET
dissimilarity_87_88 <- jaccard_dist_full["87", "88"] #FILTER VS PELLET+FILTER
dissimilarity_89_88 <- jaccard_dist_full["89", "88"] #PELLET VS PELLET+FILTER

dissimilarity_90_92 <- jaccard_dist_full["90", "92"] #FILTER VS PELLET
dissimilarity_90_91 <- jaccard_dist_full["90", "91"] #FILTER VS PELLET+FILTER
dissimilarity_92_91 <- jaccard_dist_full["92", "91"] #PELLET VS PELLET+FILTER
dissimilarity_92_93 <- jaccard_dist_full["92", "93"] #PELLET VS PELLET+MGCL2

dissimilarity_94_96 <- jaccard_dist_full["94", "96"] #FILTER VS PELLET
dissimilarity_94_95 <- jaccard_dist_full["94", "95"] #FILTER VS PELLET+FILTER
dissimilarity_96_95 <- jaccard_dist_full["96", "95"] #PELLET VS PELLET+FILTER

# Create a dataframe to nicely format the extracted dissimilarities
dissimilarity_values <- data.frame(
  Sample_Pair = c("87-89", "87-88", "89-88", "90-92", "90-91", "92-91", "92-93", "94-96", "94-95", "96-95"),
  Jaccard_Dissimilarity = c(dissimilarity_87_89, dissimilarity_87_88, dissimilarity_89_88, 
                         dissimilarity_90_92, dissimilarity_90_91, dissimilarity_92_91, 
                         dissimilarity_92_93, dissimilarity_94_96, dissimilarity_94_95, dissimilarity_96_95)
)

# Print the dissimilarity values
print(dissimilarity_values)

############################################## Diversity plots ###############################3

# Load required packages
library(ggplot2)
library(reshape2)
install.packages("lifecycle")
install.packages("dplyr")
library(dplyr)

install.packages("lifecycle", type = "binary")
install.packages("pillar", type = "binary")
install.packages("tidyselect", type = "binary")
install.packages("dplyr", type = "binary")

install.packages("rlang", type = "binary")
install.packages("vctrs", type = "binary")

# Rarefy the data to an even depth of 2400 reads/sample
ps_rarefied <- rarefy_even_depth(ps, sample.size = 2400, rngseed = 1)

# Subset the phyloseq object to keep only the OTUs table
otu_table_transposed <- t(otu_table(ps_rarefied))

# Convert the OTU table to a data frame
df_otu <- as.data.frame(otu_table_transposed)

# Add a new column with the sample IDs
df_otu$sample_id <- sample_names(ps)

# Melt the data frame to long format
df_melted <- melt(df_otu, id.vars = "sample_id")
colnames(df_melted)[2] <- "ASV"

# Join the taxonomy table to the OTU table
df_taxonomy <- as.data.frame(tax_table(ps))
df_taxonomy$ASV <- rownames(df_taxonomy) # create a new column with ASV information

df_join <- left_join(df_melted, df_taxonomy, by = "ASV")

# Aggregate the relative abundance by taxonomic level and sample
df_join_agg <- aggregate(value ~ sample_id + ASV + Kingdom + Phylum + Class + Order + Family + Genus, data = df_join, sum)

# Reorder the taxonomic levels based on abundance
df_join_agg <- df_join_agg[order(-df_join_agg$value),]


# Create stacked bar plot of relative abundance by taxonomic level and sample
ggplot(df_join_agg, aes(x = sample_id, y = value, fill = value)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Sample ID") +
  ylab("Relative Abundance") +
  ggtitle("Relative Abundance of Taxonomic Levels by Sample")

############################################################################## Genus #########

library(dplyr)
library(ggplot2)

# Group the data by Genus, then summarize by the sum of value
df_genus <- df_join_agg %>%
  group_by(Genus) %>%
  summarize(sum_value = sum(value)) %>%
  ungroup()

# Identify the top 10 Genera based on the sum_value
top_10_genera <- df_genus %>%
  arrange(desc(sum_value)) %>%
  head(10) %>%
  select(Genus)

# Group the data by sample_id and Genus, then summarize by the sum of value
df_genus <- df_join_agg %>%
  group_by(sample_id, Genus) %>%
  summarize(sum_value = sum(value)) %>%
  ungroup()

# Create a new column "Genus_grouped" to group Genera into "Other" if they are not in the top 10
df_genus <- df_genus %>%
  mutate(Genus_grouped = ifelse(Genus %in% top_10_genera$Genus, Genus, "Other"))

# Keep only the columns needed for the final plot
df_genus <- df_genus %>%
  select(sample_id, Genus_grouped, sum_value)

# Summarize the data by sample_id and Genus_grouped
df_genus_combined <- df_genus %>%
  group_by(sample_id, Genus_grouped) %>%
  summarize(sum_value = sum(sum_value)) %>%
  ungroup()

# Calculate the percentage of each Genus_grouped within each sample
df_genus_combined <- df_genus_combined %>%
  group_by(sample_id) %>%
  mutate(percent_value = sum_value / sum(sum_value) * 100) %>%
  ungroup()

# Create a varied and colorful palette of colors
my_palette <- rainbow(length(unique(df_genus_combined$Genus_grouped)))

############################## Genus Auteuil #############

#Subset metadata to keep only samples from Auteuil wastewater station
Auteuil_metadata <- Metadata[Metadata$POINT == "Auteuil", ]

# Extract sample ID from row names
Auteuil_metadata$sample_id <- rownames(Auteuil_metadata)

# Remove the prefix "FR-" from the sample IDs
Auteuil_metadata$sample_id <- gsub("FR-", "", Auteuil_metadata$sample_id)

# Select only necessary columns
Auteuil_metadata <- Auteuil_metadata[, c("sample_id", "Date", "season")]

# Merge metadata with df_genus_combined
Auteuil_merged <- merge(Auteuil_metadata, df_genus_combined, by = "sample_id")

# Create the plot with the modifications and updated palette
ggplot(Auteuil_merged, aes(x = sample_id, y = percent_value, fill = Genus_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Sample ID") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))


# Group the samples by season
Auteuil_merged$season <- factor(Auteuil_merged$season,
                                        levels = c("Winter", "Spring", "Summer", "Fall"),
                                        ordered = TRUE)
Auteuil_merged <- Auteuil_merged[order(Auteuil_merged$season), ]

# Filter out the "summer" season (just 1 sample)
Auteuil_merged_filtered <- Auteuil_merged[Auteuil_merged$season != "Summer",]

# Create the plot (each season/each sample)
ggplot(Auteuil_merged_filtered, aes(x = sample_id, y = percent_value, fill = Genus_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +  
  facet_wrap(~ season, scales = "free_x") +
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Auteuil samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


################################## Genus Lapiniere ############

# Subset metadata to keep only samples from La Piniere wastewater station
# Correct selection of metadata columns
Lapiniere_metadata <- Metadata[Metadata$POINT == "La Piniere", ]
Lapiniere_metadata$sample_id <- rownames(Lapiniere_metadata)
Lapiniere_metadata$sample_id <- gsub("FR-", "", Lapiniere_metadata$sample_id)
Lapiniere_metadata <- Lapiniere_metadata[, c("sample_id", "Date", "season")]

# Merge metadata with df_genus_combined
Lapiniere_merged <- merge(Lapiniere_metadata, df_genus_combined, by = "sample_id")

# Group the data by season
Lapiniere_merged$season <- factor(Lapiniere_merged$season,
                                  levels = c("Winter", "Spring", "Summer", "Fall"),
                                  ordered = TRUE)
Lapiniere_merged <- Lapiniere_merged[order(Lapiniere_merged$season), ]

# Filter out the "Summer" season
Lapiniere_merged_filtered <- Lapiniere_merged[Lapiniere_merged$season != "Summer",]

# Create the plot considering Genus_grouped and the seasons
ggplot(Lapiniere_merged_filtered, aes(x = sample_id, y = percent_value, fill = Genus_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +  # Use the same color palette
  facet_wrap(~ season, scales = "free_x") +
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Lapiniere samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


################################## Genus Fabreville ######

# Subset metadata to keep only samples from Fabreville wastewater station
Fabreville_metadata <- Metadata[Metadata$POINT == "Fabreville", ]

# Extract sample ID from row names
Fabreville_metadata$sample_id <- rownames(Fabreville_metadata)

# Remove the prefix "FR-" from the sample IDs
Fabreville_metadata$sample_id <- gsub("FR-", "", Fabreville_metadata$sample_id)

# Select only necessary columns
Fabreville_metadata <- Fabreville_metadata[, c("sample_id", "Date", "season")]

# Merge metadata with df_genus_combined
Fabreville_merged <- merge(Fabreville_metadata, df_genus_combined, by = "sample_id")

# Group the data by season
Fabreville_merged$season <- factor(Fabreville_merged$season,
                                   levels = c("Winter", "Spring", "Summer", "Fall"),
                                   ordered = TRUE)
Fabreville_merged <- Fabreville_merged[order(Fabreville_merged$season), ]

# Filter out the "Summer" season
Fabreville_merged_filtered <- Fabreville_merged[Fabreville_merged$season != "Summer",]

# Create the plot for each season/each sample
ggplot(Fabreville_merged_filtered, aes(x = sample_id, y = percent_value, fill = Genus_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +  # Use the same color palette
  facet_wrap(~ season, scales = "free_x") +
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Fabreville samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

##################### Genus Saisons ###########

### Subseting samples acording seasons (all WWTPs)

### Subsetting samples according to seasons (all WWTPs)

# Subset metadata to keep only samples from Fall
Fall_metadata <- Metadata[Metadata$season == "Fall", ]

# Extract sample ID from row names
Fall_metadata$sample_id <- rownames(Fall_metadata)

# Remove the prefix "FR-" from the sample IDs
Fall_metadata$sample_id <- gsub("FR-", "", Fall_metadata$sample_id)

# Select only necessary columns
Fall_metadata <- Fall_metadata[, c("sample_id", "Date", "season")]

# Merge metadata with df_genus_combined
Fall_merged <- merge(Fall_metadata, df_genus_combined, by = "sample_id")

# Create the plot for all samples from Fall
ggplot(Fall_merged, aes(x = sample_id, y = percent_value, fill = Genus_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +  # Use the same color palette
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Sample ID") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

# Subset metadata to keep only samples from Spring
Spring_metadata <- Metadata[Metadata$season == "Spring", ]

# Extract sample ID from row names
Spring_metadata$sample_id <- rownames(Spring_metadata)

# Remove the prefix "FR-" from the sample IDs
Spring_metadata$sample_id <- gsub("FR-", "", Spring_metadata$sample_id)

# Select only necessary columns
Spring_metadata <- Spring_metadata[, c("sample_id", "Date", "season")]

# Merge metadata with df_genus_combined
Spring_merged <- merge(Spring_metadata, df_genus_combined, by = "sample_id")

# Create the plot for all samples from Spring
ggplot(Spring_merged, aes(x = sample_id, y = percent_value, fill = Genus_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +  # Use the same color palette
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Sample ID") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

### Subsetting samples according to seasons (all WWTPs)

# Subset metadata to keep only samples from Winter
Winter_metadata <- Metadata[Metadata$season == "Winter", ]

# Extract sample ID from row names
Winter_metadata$sample_id <- rownames(Winter_metadata)

# Remove the prefix "FR-" from the sample IDs
Winter_metadata$sample_id <- gsub("FR-", "", Winter_metadata$sample_id)

# Select only necessary columns
Winter_metadata <- Winter_metadata[, c("sample_id", "Date", "season")]

# Merge metadata with df_genus_combined
Winter_merged <- merge(Winter_metadata, df_genus_combined, by = "sample_id")

# Create the plot for all samples from Winter
ggplot(Winter_merged, aes(x = sample_id, y = percent_value, fill = Genus_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +  # Use the same color palette
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Sample ID") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

################ Genus Filters and Pellets #######

### Subsetting samples according to concentration methods (all WWTPs)

# Subset metadata to keep only samples processed with Culot
Culot_metadata <- Metadata[Metadata$MET == "Pellet", ]

# Extract sample ID from row names
Culot_metadata$sample_id <- rownames(Culot_metadata)

# Remove the prefix "FR-" from the sample IDs
Culot_metadata$sample_id <- gsub("FR-", "", Culot_metadata$sample_id)

# Select only necessary columns
Culot_metadata <- Culot_metadata[, c("sample_id", "Date", "MET")]

# Merge metadata with df_genus_combined
Culot_merged <- merge(Culot_metadata, df_genus_combined, by = "sample_id")

# Create the plot for all samples processed with Culot
ggplot(Culot_merged, aes(x = sample_id, y = percent_value, fill = Genus_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +  # Use the same color palette
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Pellet Samples") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))


### Subsetting samples according to concentration methods (all WWTPs)

# Subset metadata to keep only samples processed with Filter
Filtre_metadata <- Metadata[Metadata$MET == "Filter", ]

# Extract sample ID from row names
Filtre_metadata$sample_id <- rownames(Filtre_metadata)

# Remove the prefix "FR-" from the sample IDs
Filtre_metadata$sample_id <- gsub("FR-", "", Filtre_metadata$sample_id)

# Select only necessary columns
Filtre_metadata <- Filtre_metadata[, c("sample_id", "Date", "MET")]

# Merge metadata with df_genus_combined
Filtre_merged <- merge(Filtre_metadata, df_genus_combined, by = "sample_id")

# Create the plot for all samples processed with Filter
ggplot(Filtre_merged, aes(x = sample_id, y = percent_value, fill = Genus_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +  # Use the same color palette
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Filter Samples") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

### Subsetting samples according to concentration methods (all WWTPs)

# Subset metadata to keep only samples processed with CulotFiltre
CulotFiltre_metadata <- Metadata[Metadata$MET == "Filter and pellet", ]

# Extract sample ID from row names
CulotFiltre_metadata$sample_id <- rownames(CulotFiltre_metadata)

# Remove the prefix "FR-" from the sample IDs
CulotFiltre_metadata$sample_id <- gsub("FR-", "", CulotFiltre_metadata$sample_id)

# Select only necessary columns
CulotFiltre_metadata <- CulotFiltre_metadata[, c("sample_id", "Date", "MET")]

# Merge metadata with df_genus_combined
CulotFiltre_merged <- merge(CulotFiltre_metadata, df_genus_combined, by = "sample_id")

# Create the plot for all samples processed with CulotFiltre
ggplot(CulotFiltre_merged, aes(x = sample_id, y = percent_value, fill = Genus_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +  # Use the same color palette
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Pellet+Filter Samples") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

###################### Genus Extraction Kit #############
### Subsetting samples according to extraction kit (all WWTPs)

# Subset metadata to keep only samples processed with PowerSoil Pro
PowerSoilPro_metadata <- Metadata[Metadata$KIT == "PowerSoil Pro", ]

# Extract sample ID from row names
PowerSoilPro_metadata$sample_id <- rownames(PowerSoilPro_metadata)

# Remove the prefix "FR-" from the sample IDs
PowerSoilPro_metadata$sample_id <- gsub("FR-", "", PowerSoilPro_metadata$sample_id)

# Select only necessary columns
PowerSoilPro_metadata <- PowerSoilPro_metadata[, c("sample_id", "Date", "KIT")]

# Merge metadata with df_genus_combined
PowerSoilPro_merged <- merge(PowerSoilPro_metadata, df_genus_combined, by = "sample_id")

# Create the plot for all samples processed with PowerSoil Pro
ggplot(PowerSoilPro_merged, aes(x = sample_id, y = percent_value, fill = Genus_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +  # Use the custom color palette
  theme_classic() +
  ylab("Relative abundance") +
  xlab("PowerSoil Pro Samples") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))


### Subsetting samples according to extraction kit (all WWTPs)

# Subset metadata to keep only samples processed with PowerLyzer
PowerLyzer_metadata <- Metadata[Metadata$KIT == "PowerLyzer", ]

# Extract sample ID from row names
PowerLyzer_metadata$sample_id <- rownames(PowerLyzer_metadata)

# Remove the prefix "FR-" from the sample IDs
PowerLyzer_metadata$sample_id <- gsub("FR-", "", PowerLyzer_metadata$sample_id)

# Select only necessary columns
PowerLyzer_metadata <- PowerLyzer_metadata[, c("sample_id", "Date", "KIT")]

# Merge metadata with df_genus_combined
PowerLyzer_merged <- merge(PowerLyzer_metadata, df_genus_combined, by = "sample_id")

# Create the plot for all samples processed with PowerLyzer
ggplot(PowerLyzer_merged, aes(x = sample_id, y = percent_value, fill = Genus_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +  # Use the same color palette
  theme_classic() +
  ylab("Relative abundance") +
  xlab("PowerLyzer Samples") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

########################### Genus Volume #######################
### Subsetting samples according to volume (all WWTPs)

# Subset metadata to keep only samples with a volume of 50
Volume50_metadata <- Metadata[Metadata$VOL == 50, ]

# Extract sample ID from row names
Volume50_metadata$sample_id <- rownames(Volume50_metadata)

# Remove the prefix "FR-" from the sample IDs
Volume50_metadata$sample_id <- gsub("FR-", "", Volume50_metadata$sample_id)

# Select only necessary columns
Volume50_metadata <- Volume50_metadata[, c("sample_id", "Date", "VOL")]

# Merge metadata with df_genus_combined
Volume50_merged <- merge(Volume50_metadata, df_genus_combined, by = "sample_id")

# Create the plot for all samples with volume of 50
ggplot(Volume50_merged, aes(x = sample_id, y = percent_value, fill = Genus_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +  # Use the same color palette
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Volume 50 Samples") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))


# Combining metadata for both volumes 100 and 125
Volume100_metadata <- Metadata[Metadata$VOL == 100, ]
Volume125_metadata <- Metadata[Metadata$VOL == 125, ]

# Prepare both subsets
Volume100_metadata$sample_id <- gsub("FR-", "", rownames(Volume100_metadata))
Volume125_metadata$sample_id <- gsub("FR-", "", rownames(Volume125_metadata))

# Combine the two datasets
Combined_metadata <- rbind(Volume100_metadata, Volume125_metadata)
Combined_metadata <- Combined_metadata[, c("sample_id", "Date", "VOL")]

# Merge combined metadata with df_genus_combined
Combined_merged <- merge(Combined_metadata, df_genus_combined, by = "sample_id")

# Create the plot for samples with volumes of 100 and 125 together
ggplot(Combined_merged, aes(x = sample_id, y = percent_value, fill = Genus_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Sample ID") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

############## Genus Centrifugation Speeed ##########################

# Subset metadata to keep only samples with a centrifugation speed of 3400
Speed3400_metadata <- Metadata[Metadata$VITESSE == 3400, ]

# Extract sample ID from row names
Speed3400_metadata$sample_id <- gsub("FR-", "", rownames(Speed3400_metadata))

# Select only necessary columns
Speed3400_metadata <- Speed3400_metadata[, c("sample_id", "Date", "VITESSE")]

# Merge metadata with df_genus_combined
Speed3400_merged <- merge(Speed3400_metadata, df_genus_combined, by = "sample_id")

# Create the plot for all samples with a centrifugation speed of 3400
ggplot(Speed3400_merged, aes(x = sample_id, y = percent_value, fill = Genus_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Centrifugation Speed 3400 Samples") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Subset metadata to keep only samples with a centrifugation speed of 15000
Speed15000_metadata <- Metadata[Metadata$VITESSE == 15000, ]

# Extract sample ID from row names
Speed15000_metadata$sample_id <- gsub("FR-", "", rownames(Speed15000_metadata))

# Select only necessary columns
Speed15000_metadata <- Speed15000_metadata[, c("sample_id", "Date", "VITESSE")]

# Merge metadata with df_genus_combined
Speed15000_merged <- merge(Speed15000_metadata, df_genus_combined, by = "sample_id")

# Create the plot for all samples with a centrifugation speed of 15000
ggplot(Speed15000_merged, aes(x = sample_id, y = percent_value, fill = Genus_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Centrifugation Speed 15000 Samples") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


##################################Genus Volume Pretest #########

# Filter metadata to keep only samples from the specific date
DateSpecific_metadata <- Metadata[Metadata$Date == "2021-12-07", ]

# Extract sample ID from row names, assuming rownames contain the sample IDs
DateSpecific_metadata$sample_id <- rownames(DateSpecific_metadata)

# Select only the necessary columns
DateSpecific_metadata <- DateSpecific_metadata[, c("sample_id", "Date")]

# Merge metadata with df_genus_combined
# Ensure that the ID column in df_genus_combined matches the sample_id
DateSpecific_merged <- merge(DateSpecific_metadata, df_genus_combined, by = "sample_id")

# Create the plot for all samples collected on the specific date
ggplot(DateSpecific_merged, aes(x = sample_id, y = percent_value, fill = Genus_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +
  theme_classic() +
  ylab("Relative Abundance") +
  xlab("Samples Collected on 2021-12-07") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

######################################## Genus Extraction Pretest##############

# Filter metadata to keep only samples from the specific date
DateSpecific_metadata <- Metadata[Metadata$Date == "2022-01-20", ]

# Extract sample ID from row names, assuming rownames contain the sample IDs
DateSpecific_metadata$sample_id <- rownames(DateSpecific_metadata)

# Select only the necessary columns
DateSpecific_metadata <- DateSpecific_metadata[, c("sample_id", "Date")]

# Merge metadata with df_genus_combined
# Ensure that the ID column in df_genus_combined matches the sample_id
DateSpecific_merged <- merge(DateSpecific_metadata, df_genus_combined, by = "sample_id")

# Create the plot for all samples collected on the specific date
ggplot(DateSpecific_merged, aes(x = sample_id, y = percent_value, fill = Genus_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +
  theme_classic() +
  ylab("Relative Abundance") +
  xlab("Samples Collected on 2022-01-20") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


############################### Genus Speed Pretest #######################3

# Filter metadata to keep only samples from the specific date
DateSpecific_metadata <- Metadata[Metadata$Date == "2022-02-17", ]

# Extract sample ID from row names, assuming rownames contain the sample IDs
DateSpecific_metadata$sample_id <- rownames(DateSpecific_metadata)

# Select only the necessary columns
DateSpecific_metadata <- DateSpecific_metadata[, c("sample_id", "Date")]

# Merge metadata with df_genus_combined
# Ensure that the ID column in df_genus_combined matches the sample_id
DateSpecific_merged <- merge(DateSpecific_metadata, df_genus_combined, by = "sample_id")

# Create the plot for all samples collected on the specific date
ggplot(DateSpecific_merged, aes(x = sample_id, y = percent_value, fill = Genus_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +
  theme_classic() +
  ylab("Relative Abundance") +
  xlab("Samples Collected on 2021-12-07") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

################################### Genus Filter and Pellet Pretest #########

# Filter metadata to keep only samples from the specific date
DateSpecific_metadata <- Metadata[Metadata$Date == "2022-04-21", ]

# Extract sample ID from row names, assuming rownames contain the sample IDs
DateSpecific_metadata$sample_id <- rownames(DateSpecific_metadata)

# Select only the necessary columns
DateSpecific_metadata <- DateSpecific_metadata[, c("sample_id", "Date")]

# Merge metadata with df_genus_combined
# Ensure that the ID column in df_genus_combined matches the sample_id
DateSpecific_merged <- merge(DateSpecific_metadata, df_genus_combined, by = "sample_id")

# Create the plot for all samples collected on the specific date
ggplot(DateSpecific_merged, aes(x = sample_id, y = percent_value, fill = Genus_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +
  theme_classic() +
  ylab("Relative Abundance") +
  xlab("Samples Collected on 2021-12-07") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

###################################################Analysis by Class###########

#######10 dominant Class

# Group the data by Class, then summarize by the sum of value
df_class <- df_join_agg %>%
  group_by(Class) %>%
  summarize(sum_value = sum(value)) %>%
  ungroup()

# Identify the top 10 Classes based on the sum_value
top_10_classes <- df_class %>%
  arrange(desc(sum_value)) %>%
  head(10) %>%
  select(Class)

# Group the data by sample_id and Class, then summarize by the sum of value
df_class <- df_join_agg %>%
  group_by(sample_id, Class) %>%
  summarize(sum_value = sum(value)) %>%
  ungroup()

# Create a new column "Class_grouped" to group Classes into "Other" if they are not in the top 10
df_class <- df_class %>%
  mutate(Class_grouped = ifelse(Class %in% top_10_classes$Class, Class, "Other"))

# Keep only the columns needed for the final plot
df_class <- df_class %>%
  select(sample_id, Class_grouped, sum_value)

# Summarize the data by sample_id and Class_grouped
df_class_combined <- df_class %>%
  group_by(sample_id, Class_grouped) %>%
  summarize(sum_value = sum(sum_value)) %>%
  ungroup()

# Calculate the percentage of each Class_grouped within each sample
df_class_combined <- df_class_combined %>%
  group_by(sample_id) %>%
  mutate(percent_value = sum_value / sum(sum_value) * 100) %>%
  ungroup()

# Create a varied and colorful palette of colors
my_palette <- rainbow(length(unique(df_class_combined$Class_grouped)))


########################### Class Auteuil ###############


# Merge metadata with df_class_combined
Auteuil_merged <- merge(Auteuil_metadata, df_class_combined, by = "sample_id")

# Create the plot with modifications for class
ggplot(Auteuil_merged, aes(x = sample_id, y = percent_value, fill = Class_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Sample ID") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))


# Group the samples by season
Auteuil_merged$season <- factor(Auteuil_merged$season,
                                levels = c("Winter", "Spring", "Summer", "Fall"),
                                ordered = TRUE)
Auteuil_merged <- Auteuil_merged[order(Auteuil_merged$season), ]

# Filter out the "summer" season (just 1 sample)
Auteuil_merged_filtered <- Auteuil_merged[Auteuil_merged$season != "Summer",]

# Create the plot (each season/each sample)
ggplot(Auteuil_merged_filtered, aes(x = sample_id, y = percent_value, fill = Class_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +  # Usar a mesma paleta de cores do primeiro script
  facet_wrap(~ season, scales = "free_x") +
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Auteuil samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


################################## Class Lapiniere ###############3

# Subset metadata to keep only samples from La Piniere wastewater station
# Correct selection of metadata columns
Lapiniere_metadata <- Metadata[Metadata$POINT == "La Piniere", ]
Lapiniere_metadata$sample_id <- rownames(Lapiniere_metadata)
Lapiniere_metadata$sample_id <- gsub("FR-", "", Lapiniere_metadata$sample_id)
Lapiniere_metadata <- Lapiniere_metadata[, c("sample_id", "Date", "season")]

# Merge metadata with df_genus_combined
Lapiniere_merged <- merge(Lapiniere_metadata, df_class_combined, by = "sample_id")

# Group the data by season
Lapiniere_merged$season <- factor(Lapiniere_merged$season,
                                  levels = c("Winter", "Spring", "Summer", "Fall"),
                                  ordered = TRUE)
Lapiniere_merged <- Lapiniere_merged[order(Lapiniere_merged$season), ]

# Filter out the "Summer" season
Lapiniere_merged_filtered <- Lapiniere_merged[Lapiniere_merged$season != "Summer",]

# Create the plot considering Genus_grouped and the seasons
ggplot(Lapiniere_merged_filtered, aes(x = sample_id, y = percent_value, fill = Class_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +  # Use the same color palette
  facet_wrap(~ season, scales = "free_x") +
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Lapiniere samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


################################## Class Fabreville ############

# Subset metadata to keep only samples from Fabreville wastewater station
Fabreville_metadata <- Metadata[Metadata$POINT == "Fabreville", ]

# Extract sample ID from row names
Fabreville_metadata$sample_id <- rownames(Fabreville_metadata)

# Remove the prefix "FR-" from the sample IDs
Fabreville_metadata$sample_id <- gsub("FR-", "", Fabreville_metadata$sample_id)

# Select only necessary columns
Fabreville_metadata <- Fabreville_metadata[, c("sample_id", "Date", "season")]

# Merge metadata with df_genus_combined
Fabreville_merged <- merge(Fabreville_metadata, df_class_combined, by = "sample_id")

# Group the data by season
Fabreville_merged$season <- factor(Fabreville_merged$season,
                                   levels = c("Winter", "Spring", "Summer", "Fall"),
                                   ordered = TRUE)
Fabreville_merged <- Fabreville_merged[order(Fabreville_merged$season), ]

# Filter out the "Summer" season
Fabreville_merged_filtered <- Fabreville_merged[Fabreville_merged$season != "Summer",]

# Create the plot for each season/each sample
ggplot(Fabreville_merged_filtered, aes(x = sample_id, y = percent_value, fill = Class_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +  # Use the same color palette
  facet_wrap(~ season, scales = "free_x") +
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Fabreville samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

##################### Class saisons ###########

### Subseting samples acording seasons (all WWTPs)

### Subsetting samples according to seasons (all WWTPs)

# Subset metadata to keep only samples from Fall
Fall_metadata <- Metadata[Metadata$season == "Fall", ]

# Extract sample ID from row names
Fall_metadata$sample_id <- rownames(Fall_metadata)

# Remove the prefix "FR-" from the sample IDs
Fall_metadata$sample_id <- gsub("FR-", "", Fall_metadata$sample_id)

# Select only necessary columns
Fall_metadata <- Fall_metadata[, c("sample_id", "Date", "season")]

# Merge metadata with df_genus_combined
Fall_merged <- merge(Fall_metadata, df_class_combined, by = "sample_id")

# Create the plot for all samples from Fall
ggplot(Fall_merged, aes(x = sample_id, y = percent_value, fill = Class_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +  # Use the same color palette
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Sample ID") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

# Subset metadata to keep only samples from Spring
Spring_metadata <- Metadata[Metadata$season == "Spring", ]

# Extract sample ID from row names
Spring_metadata$sample_id <- rownames(Spring_metadata)

# Remove the prefix "FR-" from the sample IDs
Spring_metadata$sample_id <- gsub("FR-", "", Spring_metadata$sample_id)

# Select only necessary columns
Spring_metadata <- Spring_metadata[, c("sample_id", "Date", "season")]

# Merge metadata with df_genus_combined
Spring_merged <- merge(Spring_metadata, df_class_combined, by = "sample_id")

# Create the plot for all samples from Spring
ggplot(Spring_merged, aes(x = sample_id, y = percent_value, fill = Class_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +  # Use the same color palette
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Sample ID") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

### Subsetting samples according to seasons (all WWTPs)

# Subset metadata to keep only samples from Winter
Winter_metadata <- Metadata[Metadata$season == "Winter", ]

# Extract sample ID from row names
Winter_metadata$sample_id <- rownames(Winter_metadata)

# Remove the prefix "FR-" from the sample IDs
Winter_metadata$sample_id <- gsub("FR-", "", Winter_metadata$sample_id)

# Select only necessary columns
Winter_metadata <- Winter_metadata[, c("sample_id", "Date", "season")]

# Merge metadata with df_genus_combined
Winter_merged <- merge(Winter_metadata, df_class_combined, by = "sample_id")

# Create the plot for all samples from Winter
ggplot(Winter_merged, aes(x = sample_id, y = percent_value, fill = Class_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +  # Use the same color palette
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Sample ID") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))


#################### Class Filters and Pellets ################


### Subsetting samples according to concentration methods (all WWTPs)

# Subset metadata to keep only samples processed with Culot
Culot_metadata <- Metadata[Metadata$MET == "Pellet", ]

# Extract sample ID from row names
Culot_metadata$sample_id <- rownames(Culot_metadata)

# Remove the prefix "FR-" from the sample IDs
Culot_metadata$sample_id <- gsub("FR-", "", Culot_metadata$sample_id)

# Select only necessary columns
Culot_metadata <- Culot_metadata[, c("sample_id", "Date", "MET")]

# Merge metadata with df_genus_combined
Culot_merged <- merge(Culot_metadata, df_class_combined, by = "sample_id")

# Create the plot for all samples processed with Culot
ggplot(Culot_merged, aes(x = sample_id, y = percent_value, fill = Class_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +  # Use the same color palette
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Pellet Samples") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))


### Subsetting samples according to concentration methods (all WWTPs)

# Subset metadata to keep only samples processed with Filter
Filtre_metadata <- Metadata[Metadata$MET == "Filter", ]

# Extract sample ID from row names
Filtre_metadata$sample_id <- rownames(Filtre_metadata)

# Remove the prefix "FR-" from the sample IDs
Filtre_metadata$sample_id <- gsub("FR-", "", Filtre_metadata$sample_id)

# Select only necessary columns
Filtre_metadata <- Filtre_metadata[, c("sample_id", "Date", "MET")]

# Merge metadata with df_genus_combined
Filtre_merged <- merge(Filtre_metadata, df_class_combined, by = "sample_id")

# Create the plot for all samples processed with Filter
ggplot(Filtre_merged, aes(x = sample_id, y = percent_value, fill = Class_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +  # Use the same color palette
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Filter Samples") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

### Subsetting samples according to concentration methods (all WWTPs)

# Subset metadata to keep only samples processed with CulotFiltre
CulotFiltre_metadata <- Metadata[Metadata$MET == "Filter and pellet", ]

# Extract sample ID from row names
CulotFiltre_metadata$sample_id <- rownames(CulotFiltre_metadata)

# Remove the prefix "FR-" from the sample IDs
CulotFiltre_metadata$sample_id <- gsub("FR-", "", CulotFiltre_metadata$sample_id)

# Select only necessary columns
CulotFiltre_metadata <- CulotFiltre_metadata[, c("sample_id", "Date", "MET")]

# Merge metadata with df_genus_combined
CulotFiltre_merged <- merge(CulotFiltre_metadata, df_class_combined, by = "sample_id")

# Create the plot for all samples processed with CulotFiltre
ggplot(CulotFiltre_merged, aes(x = sample_id, y = percent_value, fill = Class_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +  # Use the same color palette
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Pellet+Filter Samples") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))


###################### Class extraction kits ############


### Subsetting samples according to extraction kit (all WWTPs)

# Subset metadata to keep only samples processed with PowerSoil Pro
PowerSoilPro_metadata <- Metadata[Metadata$KIT == "PowerSoil Pro", ]

# Extract sample ID from row names
PowerSoilPro_metadata$sample_id <- rownames(PowerSoilPro_metadata)

# Remove the prefix "FR-" from the sample IDs
PowerSoilPro_metadata$sample_id <- gsub("FR-", "", PowerSoilPro_metadata$sample_id)

# Select only necessary columns
PowerSoilPro_metadata <- PowerSoilPro_metadata[, c("sample_id", "Date", "KIT")]

# Merge metadata with df_genus_combined
PowerSoilPro_merged <- merge(PowerSoilPro_metadata, df_class_combined, by = "sample_id")

# Create the plot for all samples processed with PowerSoil Pro
ggplot(PowerSoilPro_merged, aes(x = sample_id, y = percent_value, fill = Class_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +  # Use the custom color palette
  theme_classic() +
  ylab("Relative abundance") +
  xlab("PowerSoil Pro Samples") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))


### Subsetting samples according to extraction kit (all WWTPs)

# Subset metadata to keep only samples processed with PowerLyzer
PowerLyzer_metadata <- Metadata[Metadata$KIT == "PowerLyzer", ]

# Extract sample ID from row names
PowerLyzer_metadata$sample_id <- rownames(PowerLyzer_metadata)

# Remove the prefix "FR-" from the sample IDs
PowerLyzer_metadata$sample_id <- gsub("FR-", "", PowerLyzer_metadata$sample_id)

# Select only necessary columns
PowerLyzer_metadata <- PowerLyzer_metadata[, c("sample_id", "Date", "KIT")]

# Merge metadata with df_genus_combined
PowerLyzer_merged <- merge(PowerLyzer_metadata, df_class_combined, by = "sample_id")

# Create the plot for all samples processed with PowerLyzer
ggplot(PowerLyzer_merged, aes(x = sample_id, y = percent_value, fill = Class_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +  # Use the same color palette
  theme_classic() +
  ylab("Relative abundance") +
  xlab("PowerLyzer Samples") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

################# Class Volume ######################3


### Subsetting samples according to volume (all WWTPs)

# Subset metadata to keep only samples with a volume of 50
Volume50_metadata <- Metadata[Metadata$VOL == 50, ]

# Extract sample ID from row names
Volume50_metadata$sample_id <- rownames(Volume50_metadata)

# Remove the prefix "FR-" from the sample IDs
Volume50_metadata$sample_id <- gsub("FR-", "", Volume50_metadata$sample_id)

# Select only necessary columns
Volume50_metadata <- Volume50_metadata[, c("sample_id", "Date", "VOL")]

# Merge metadata with df_genus_combined
Volume50_merged <- merge(Volume50_metadata, df_class_combined, by = "sample_id")

# Create the plot for all samples with volume of 50
ggplot(Volume50_merged, aes(x = sample_id, y = percent_value, fill = Class_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +  # Use the same color palette
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Volume 50 Samples") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))


# Combining metadata for both volumes 100 and 125
Volume100_metadata <- Metadata[Metadata$VOL == 100, ]
Volume125_metadata <- Metadata[Metadata$VOL == 125, ]

# Prepare both subsets
Volume100_metadata$sample_id <- gsub("FR-", "", rownames(Volume100_metadata))
Volume125_metadata$sample_id <- gsub("FR-", "", rownames(Volume125_metadata))

# Combine the two datasets
Combined_metadata <- rbind(Volume100_metadata, Volume125_metadata)
Combined_metadata <- Combined_metadata[, c("sample_id", "Date", "VOL")]

# Merge combined metadata with df_genus_combined
Combined_merged <- merge(Combined_metadata, df_class_combined, by = "sample_id")

# Create the plot for samples with volumes of 100 and 125 together
ggplot(Combined_merged, aes(x = sample_id, y = percent_value, fill = Class_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Sample ID") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

######################### Class Centrifugation Speed #####################

# Subset metadata to keep only samples with a centrifugation speed of 3400
Speed3400_metadata <- Metadata[Metadata$VITESSE == 3400, ]

# Extract sample ID from row names
Speed3400_metadata$sample_id <- gsub("FR-", "", rownames(Speed3400_metadata))

# Select only necessary columns
Speed3400_metadata <- Speed3400_metadata[, c("sample_id", "Date", "VITESSE")]

# Merge metadata with df_genus_combined
Speed3400_merged <- merge(Speed3400_metadata, df_class_combined, by = "sample_id")

# Create the plot for all samples with a centrifugation speed of 3400
ggplot(Speed3400_merged, aes(x = sample_id, y = percent_value, fill = Class_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Centrifugation Speed 3400 Samples") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Subset metadata to keep only samples with a centrifugation speed of 15000
Speed15000_metadata <- Metadata[Metadata$VITESSE == 15000, ]

# Extract sample ID from row names
Speed15000_metadata$sample_id <- gsub("FR-", "", rownames(Speed15000_metadata))

# Select only necessary columns
Speed15000_metadata <- Speed15000_metadata[, c("sample_id", "Date", "VITESSE")]

# Merge metadata with df_genus_combined
Speed15000_merged <- merge(Speed15000_metadata, df_class_combined, by = "sample_id")

# Create the plot for all samples with a centrifugation speed of 15000
ggplot(Speed15000_merged, aes(x = sample_id, y = percent_value, fill = Class_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Centrifugation Speed 15000 Samples") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


################################## Class Volume Pretest ######################

# Filter metadata to keep only samples from the specific date
DateSpecific_metadata <- Metadata[Metadata$Date == "2021-12-07", ]

# Extract sample ID from row names, assuming rownames contain the sample IDs
DateSpecific_metadata$sample_id <- rownames(DateSpecific_metadata)

# Select only the necessary columns
DateSpecific_metadata <- DateSpecific_metadata[, c("sample_id", "Date")]

# Merge metadata with df_genus_combined
# Ensure that the ID column in df_genus_combined matches the sample_id
DateSpecific_merged <- merge(DateSpecific_metadata, df_class_combined, by = "sample_id")

# Create the plot for all samples collected on the specific date
ggplot(DateSpecific_merged, aes(x = sample_id, y = percent_value, fill = Class_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +
  theme_classic() +
  ylab("Relative Abundance") +
  xlab("Samples Collected on 2021-12-07") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

######################################## Class Extraction Pretest #####################

# Filter metadata to keep only samples from the specific date
DateSpecific_metadata <- Metadata[Metadata$Date == "2022-01-20", ]

# Extract sample ID from row names, assuming rownames contain the sample IDs
DateSpecific_metadata$sample_id <- rownames(DateSpecific_metadata)

# Select only the necessary columns
DateSpecific_metadata <- DateSpecific_metadata[, c("sample_id", "Date")]

# Merge metadata with df_genus_combined
# Ensure that the ID column in df_genus_combined matches the sample_id
DateSpecific_merged <- merge(DateSpecific_metadata, df_class_combined, by = "sample_id")

# Create the plot for all samples collected on the specific date
ggplot(DateSpecific_merged, aes(x = sample_id, y = percent_value, fill = Class_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +
  theme_classic() +
  ylab("Relative Abundance") +
  xlab("Samples Collected on 2022-01-20") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


############################### Class Speed Pretest ###############################

# Filter metadata to keep only samples from the specific date
DateSpecific_metadata <- Metadata[Metadata$Date == "2022-02-17", ]

# Extract sample ID from row names, assuming rownames contain the sample IDs
DateSpecific_metadata$sample_id <- rownames(DateSpecific_metadata)

# Select only the necessary columns
DateSpecific_metadata <- DateSpecific_metadata[, c("sample_id", "Date")]

# Merge metadata with df_genus_combined
# Ensure that the ID column in df_genus_combined matches the sample_id
DateSpecific_merged <- merge(DateSpecific_metadata, df_class_combined, by = "sample_id")

# Create the plot for all samples collected on the specific date
ggplot(DateSpecific_merged, aes(x = sample_id, y = percent_value, fill = Class_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +
  theme_classic() +
  ylab("Relative Abundance") +
  xlab("Samples Collected on 2021-12-07") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

################################### Class Filters and Pellets Pretest ################

# Filter metadata to keep only samples from the specific date
DateSpecific_metadata <- Metadata[Metadata$Date == "2022-04-21", ]

# Extract sample ID from row names, assuming rownames contain the sample IDs
DateSpecific_metadata$sample_id <- rownames(DateSpecific_metadata)

# Select only the necessary columns
DateSpecific_metadata <- DateSpecific_metadata[, c("sample_id", "Date")]

# Merge metadata with df_genus_combined
# Ensure that the ID column in df_genus_combined matches the sample_id
DateSpecific_merged <- merge(DateSpecific_metadata, df_class_combined, by = "sample_id")

# Create the plot for all samples collected on the specific date
ggplot(DateSpecific_merged, aes(x = sample_id, y = percent_value, fill = Class_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +
  theme_classic() +
  ylab("Relative Abundance") +
  xlab("Samples Collected on 2021-12-07") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

############################################ Family #########################3
library(dplyr)
library(ggplot2)

# Group the data by Family, then summarize by the sum of value
df_family <- df_join_agg %>%
  group_by(Family) %>%
  summarize(sum_value = sum(value)) %>%
  ungroup()

# Identify the top 10 Families based on the sum_value
top_10_families <- df_family %>%
  arrange(desc(sum_value)) %>%
  head(10) %>%
  select(Family)

# Group the data by sample_id and Family, then summarize by the sum of value
df_family <- df_join_agg %>%
  group_by(sample_id, Family) %>%
  summarize(sum_value = sum(value)) %>%
  ungroup()

# Create a new column "Family_grouped" to group Families into "Other" if they are not in the top 10
df_family <- df_family %>%
  mutate(Family_grouped = ifelse(Family %in% top_10_families$Family, Family, "Other"))

# Keep only the columns needed for the final plot
df_family <- df_family %>%
  select(sample_id, Family_grouped, sum_value)

# Summarize the data by sample_id and Family_grouped
df_family_combined <- df_family %>%
  group_by(sample_id, Family_grouped) %>%
  summarize(sum_value = sum(sum_value)) %>%
  ungroup()

# Calculate the percentage of each Family_grouped within each sample
df_family_combined <- df_family_combined %>%
  group_by(sample_id) %>%
  mutate(percent_value = sum_value / sum(sum_value) * 100) %>%
  ungroup()

# Create a varied and colorful palette of colors
my_palette <- rainbow(length(unique(df_family_combined$Family_grouped)))



########################### Family Auteuil ###############


# Merge metadata with df_family_combined
Auteuil_merged <- merge(Auteuil_metadata, df_family_combined, by = "sample_id")

# Create the plot with modifications for class
ggplot(Auteuil_merged, aes(x = sample_id, y = percent_value, fill = Family_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Sample ID") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))


# Group the samples by season
Auteuil_merged$season <- factor(Auteuil_merged$season,
                                levels = c("Winter", "Spring", "Summer", "Fall"),
                                ordered = TRUE)
Auteuil_merged <- Auteuil_merged[order(Auteuil_merged$season), ]

# Filter out the "summer" season (just 1 sample)
Auteuil_merged_filtered <- Auteuil_merged[Auteuil_merged$season != "Summer",]

# Create the plot (each season/each sample)
ggplot(Auteuil_merged_filtered, aes(x = sample_id, y = percent_value, fill = Family_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +  # Usar a mesma paleta de cores do primeiro script
  facet_wrap(~ season, scales = "free_x") +
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Auteuil samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


################################## Family Lapiniere ###############3

# Subset metadata to keep only samples from La Piniere wastewater station
# Correct selection of metadata columns
Lapiniere_metadata <- Metadata[Metadata$POINT == "La Piniere", ]
Lapiniere_metadata$sample_id <- rownames(Lapiniere_metadata)
Lapiniere_metadata$sample_id <- gsub("FR-", "", Lapiniere_metadata$sample_id)
Lapiniere_metadata <- Lapiniere_metadata[, c("sample_id", "Date", "season")]

# Merge metadata with df_family_combined
Lapiniere_merged <- merge(Lapiniere_metadata, df_family_combined, by = "sample_id")

# Group the data by season
Lapiniere_merged$season <- factor(Lapiniere_merged$season,
                                  levels = c("Winter", "Spring", "Summer", "Fall"),
                                  ordered = TRUE)
Lapiniere_merged <- Lapiniere_merged[order(Lapiniere_merged$season), ]

# Filter out the "Summer" season
Lapiniere_merged_filtered <- Lapiniere_merged[Lapiniere_merged$season != "Summer",]

# Create the plot considering Family_grouped and the seasons
ggplot(Lapiniere_merged_filtered, aes(x = sample_id, y = percent_value, fill = Family_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +  # Use the same color palette
  facet_wrap(~ season, scales = "free_x") +
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Lapiniere samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


################################## Family Fabreville ############

# Subset metadata to keep only samples from Fabreville wastewater station
Fabreville_metadata <- Metadata[Metadata$POINT == "Fabreville", ]

# Extract sample ID from row names
Fabreville_metadata$sample_id <- rownames(Fabreville_metadata)

# Remove the prefix "FR-" from the sample IDs
Fabreville_metadata$sample_id <- gsub("FR-", "", Fabreville_metadata$sample_id)

# Select only necessary columns
Fabreville_metadata <- Fabreville_metadata[, c("sample_id", "Date", "season")]

# Merge metadata with df_family_combined
Fabreville_merged <- merge(Fabreville_metadata, df_family_combined, by = "sample_id")

# Group the data by season
Fabreville_merged$season <- factor(Fabreville_merged$season,
                                   levels = c("Winter", "Spring", "Summer", "Fall"),
                                   ordered = TRUE)
Fabreville_merged <- Fabreville_merged[order(Fabreville_merged$season), ]

# Filter out the "Summer" season
Fabreville_merged_filtered <- Fabreville_merged[Fabreville_merged$season != "Summer",]

# Create the plot for each season/each sample
ggplot(Fabreville_merged_filtered, aes(x = sample_id, y = percent_value, fill = Family_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +  # Use the same color palette
  facet_wrap(~ season, scales = "free_x") +
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Fabreville samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

##################### Family saisons ###########

### Subseting samples acording seasons (all WWTPs)

### Subsetting samples according to seasons (all WWTPs)

# Subset metadata to keep only samples from Fall
Fall_metadata <- Metadata[Metadata$season == "Fall", ]

# Extract sample ID from row names
Fall_metadata$sample_id <- rownames(Fall_metadata)

# Remove the prefix "FR-" from the sample IDs
Fall_metadata$sample_id <- gsub("FR-", "", Fall_metadata$sample_id)

# Select only necessary columns
Fall_metadata <- Fall_metadata[, c("sample_id", "Date", "season")]

# Merge metadata with df_genus_combined
Fall_merged <- merge(Fall_metadata, df_family_combined, by = "sample_id")

# Create the plot for all samples from Fall
ggplot(Fall_merged, aes(x = sample_id, y = percent_value, fill = Family_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +  # Use the same color palette
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Sample ID") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

# Subset metadata to keep only samples from Spring
Spring_metadata <- Metadata[Metadata$season == "Spring", ]

# Extract sample ID from row names
Spring_metadata$sample_id <- rownames(Spring_metadata)

# Remove the prefix "FR-" from the sample IDs
Spring_metadata$sample_id <- gsub("FR-", "", Spring_metadata$sample_id)

# Select only necessary columns
Spring_metadata <- Spring_metadata[, c("sample_id", "Date", "season")]

# Merge metadata with df_family_combined
Spring_merged <- merge(Spring_metadata, df_family_combined, by = "sample_id")

# Create the plot for all samples from Spring
ggplot(Spring_merged, aes(x = sample_id, y = percent_value, fill = Family_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +  # Use the same color palette
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Sample ID") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

### Subsetting samples according to seasons (all WWTPs)

# Subset metadata to keep only samples from Winter
Winter_metadata <- Metadata[Metadata$season == "Winter", ]

# Extract sample ID from row names
Winter_metadata$sample_id <- rownames(Winter_metadata)

# Remove the prefix "FR-" from the sample IDs
Winter_metadata$sample_id <- gsub("FR-", "", Winter_metadata$sample_id)

# Select only necessary columns
Winter_metadata <- Winter_metadata[, c("sample_id", "Date", "season")]

# Merge metadata with df_genus_combined
Winter_merged <- merge(Winter_metadata, df_family_combined, by = "sample_id")

# Create the plot for all samples from Winter
ggplot(Winter_merged, aes(x = sample_id, y = percent_value, fill = Family_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +  # Use the same color palette
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Sample ID") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))


#################### Family Filters and Pellets ################


### Subsetting samples according to concentration methods (all WWTPs)

# Subset metadata to keep only samples processed with Culot
Culot_metadata <- Metadata[Metadata$MET == "Pellet", ]

# Extract sample ID from row names
Culot_metadata$sample_id <- rownames(Culot_metadata)

# Remove the prefix "FR-" from the sample IDs
Culot_metadata$sample_id <- gsub("FR-", "", Culot_metadata$sample_id)

# Select only necessary columns
Culot_metadata <- Culot_metadata[, c("sample_id", "Date", "MET")]

# Merge metadata with df_family_combined
Culot_merged <- merge(Culot_metadata, df_family_combined, by = "sample_id")

# Create the plot for all samples processed with Culot
ggplot(Culot_merged, aes(x = sample_id, y = percent_value, fill = Family_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +  # Use the same color palette
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Pellet Samples") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))


### Subsetting samples according to concentration methods (all WWTPs)

# Subset metadata to keep only samples processed with Filter
Filtre_metadata <- Metadata[Metadata$MET == "Filter", ]

# Extract sample ID from row names
Filtre_metadata$sample_id <- rownames(Filtre_metadata)

# Remove the prefix "FR-" from the sample IDs
Filtre_metadata$sample_id <- gsub("FR-", "", Filtre_metadata$sample_id)

# Select only necessary columns
Filtre_metadata <- Filtre_metadata[, c("sample_id", "Date", "MET")]

# Merge metadata with df_family_combined
Filtre_merged <- merge(Filtre_metadata, df_family_combined, by = "sample_id")

# Create the plot for all samples processed with Filter
ggplot(Filtre_merged, aes(x = sample_id, y = percent_value, fill = Family_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +  # Use the same color palette
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Filter Samples") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

### Subsetting samples according to concentration methods (all WWTPs)

# Subset metadata to keep only samples processed with CulotFiltre
CulotFiltre_metadata <- Metadata[Metadata$MET == "Filter and pellet", ]

# Extract sample ID from row names
CulotFiltre_metadata$sample_id <- rownames(CulotFiltre_metadata)

# Remove the prefix "FR-" from the sample IDs
CulotFiltre_metadata$sample_id <- gsub("FR-", "", CulotFiltre_metadata$sample_id)

# Select only necessary columns
CulotFiltre_metadata <- CulotFiltre_metadata[, c("sample_id", "Date", "MET")]

# Merge metadata with df_family_combined
CulotFiltre_merged <- merge(CulotFiltre_metadata, df_family_combined, by = "sample_id")

# Create the plot for all samples processed with CulotFiltre
ggplot(CulotFiltre_merged, aes(x = sample_id, y = percent_value, fill = Family_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +  # Use the same color palette
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Pellet+Filter Samples") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))


###################### Family extraction kits ############


### Subsetting samples according to extraction kit (all WWTPs)

# Subset metadata to keep only samples processed with PowerSoil Pro
PowerSoilPro_metadata <- Metadata[Metadata$KIT == "PowerSoil Pro", ]

# Extract sample ID from row names
PowerSoilPro_metadata$sample_id <- rownames(PowerSoilPro_metadata)

# Remove the prefix "FR-" from the sample IDs
PowerSoilPro_metadata$sample_id <- gsub("FR-", "", PowerSoilPro_metadata$sample_id)

# Select only necessary columns
PowerSoilPro_metadata <- PowerSoilPro_metadata[, c("sample_id", "Date", "KIT")]

# Merge metadata with df_family_combined
PowerSoilPro_merged <- merge(PowerSoilPro_metadata, df_family_combined, by = "sample_id")

# Create the plot for all samples processed with PowerSoil Pro
ggplot(PowerSoilPro_merged, aes(x = sample_id, y = percent_value, fill = Family_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +  # Use the custom color palette
  theme_classic() +
  ylab("Relative abundance") +
  xlab("PowerSoil Pro Samples") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))


### Subsetting samples according to extraction kit (all WWTPs)

# Subset metadata to keep only samples processed with PowerLyzer
PowerLyzer_metadata <- Metadata[Metadata$KIT == "PowerLyzer", ]

# Extract sample ID from row names
PowerLyzer_metadata$sample_id <- rownames(PowerLyzer_metadata)

# Remove the prefix "FR-" from the sample IDs
PowerLyzer_metadata$sample_id <- gsub("FR-", "", PowerLyzer_metadata$sample_id)

# Select only necessary columns
PowerLyzer_metadata <- PowerLyzer_metadata[, c("sample_id", "Date", "KIT")]

# Merge metadata with df_family_combined
PowerLyzer_merged <- merge(PowerLyzer_metadata, df_family_combined, by = "sample_id")

# Create the plot for all samples processed with PowerLyzer
ggplot(PowerLyzer_merged, aes(x = sample_id, y = percent_value, fill = Family_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +  # Use the same color palette
  theme_classic() +
  ylab("Relative abundance") +
  xlab("PowerLyzer Samples") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

################# Family Volume ######################


### Subsetting samples according to volume (all WWTPs)

# Subset metadata to keep only samples with a volume of 50
Volume50_metadata <- Metadata[Metadata$VOL == 50, ]

# Extract sample ID from row names
Volume50_metadata$sample_id <- rownames(Volume50_metadata)

# Remove the prefix "FR-" from the sample IDs
Volume50_metadata$sample_id <- gsub("FR-", "", Volume50_metadata$sample_id)

# Select only necessary columns
Volume50_metadata <- Volume50_metadata[, c("sample_id", "Date", "VOL")]

# Merge metadata with df_family_combined
Volume50_merged <- merge(Volume50_metadata, df_family_combined, by = "sample_id")

# Create the plot for all samples with volume of 50
ggplot(Volume50_merged, aes(x = sample_id, y = percent_value, fill = Family_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +  # Use the same color palette
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Volume 50 Samples") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))


# Combining metadata for both volumes 100 and 125
Volume100_metadata <- Metadata[Metadata$VOL == 100, ]
Volume125_metadata <- Metadata[Metadata$VOL == 125, ]

# Prepare both subsets
Volume100_metadata$sample_id <- gsub("FR-", "", rownames(Volume100_metadata))
Volume125_metadata$sample_id <- gsub("FR-", "", rownames(Volume125_metadata))

# Combine the two datasets
Combined_metadata <- rbind(Volume100_metadata, Volume125_metadata)
Combined_metadata <- Combined_metadata[, c("sample_id", "Date", "VOL")]

# Merge combined metadata with df_family_combined
Combined_merged <- merge(Combined_metadata, df_family_combined, by = "sample_id")

# Create the plot for samples with volumes of 100 and 125 together
ggplot(Combined_merged, aes(x = sample_id, y = percent_value, fill = Family_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Sample ID") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

######################### Family Centrifugation Speed #####################

# Subset metadata to keep only samples with a centrifugation speed of 3400
Speed3400_metadata <- Metadata[Metadata$VITESSE == 3400, ]

# Extract sample ID from row names
Speed3400_metadata$sample_id <- gsub("FR-", "", rownames(Speed3400_metadata))

# Select only necessary columns
Speed3400_metadata <- Speed3400_metadata[, c("sample_id", "Date", "VITESSE")]

# Merge metadata with df_family_combined
Speed3400_merged <- merge(Speed3400_metadata, df_family_combined, by = "sample_id")

# Create the plot for all samples with a centrifugation speed of 3400
ggplot(Speed3400_merged, aes(x = sample_id, y = percent_value, fill = Family_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Centrifugation Speed 3400 Samples") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Subset metadata to keep only samples with a centrifugation speed of 15000
Speed15000_metadata <- Metadata[Metadata$VITESSE == 15000, ]

# Extract sample ID from row names
Speed15000_metadata$sample_id <- gsub("FR-", "", rownames(Speed15000_metadata))

# Select only necessary columns
Speed15000_metadata <- Speed15000_metadata[, c("sample_id", "Date", "VITESSE")]

# Merge metadata with df_family_combined
Speed15000_merged <- merge(Speed15000_metadata, df_family_combined, by = "sample_id")

# Create the plot for all samples with a centrifugation speed of 15000
ggplot(Speed15000_merged, aes(x = sample_id, y = percent_value, fill = Family_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +
  theme_classic() +
  ylab("Relative abundance") +
  xlab("Centrifugation Speed 15000 Samples") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


################################## Family Volume Pretest ######################

# Filter metadata to keep only samples from the specific date
DateSpecific_metadata <- Metadata[Metadata$Date == "2021-12-07", ]

# Extract sample ID from row names, assuming rownames contain the sample IDs
DateSpecific_metadata$sample_id <- rownames(DateSpecific_metadata)

# Select only the necessary columns
DateSpecific_metadata <- DateSpecific_metadata[, c("sample_id", "Date")]

# Merge metadata with df_family_combined
# Ensure that the ID column in df_genus_combined matches the sample_id
DateSpecific_merged <- merge(DateSpecific_metadata, df_family_combined, by = "sample_id")

# Create the plot for all samples collected on the specific date
ggplot(DateSpecific_merged, aes(x = sample_id, y = percent_value, fill = Family_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +
  theme_classic() +
  ylab("Relative Abundance") +
  xlab("Samples Collected on 2021-12-07") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

######################################## Family Extraction Pretest #####################

# Filter metadata to keep only samples from the specific date
DateSpecific_metadata <- Metadata[Metadata$Date == "2022-01-20", ]

# Extract sample ID from row names, assuming rownames contain the sample IDs
DateSpecific_metadata$sample_id <- rownames(DateSpecific_metadata)

# Select only the necessary columns
DateSpecific_metadata <- DateSpecific_metadata[, c("sample_id", "Date")]

# Merge metadata with df_genus_combined
# Ensure that the ID column in df_family_combined matches the sample_id
DateSpecific_merged <- merge(DateSpecific_metadata, df_family_combined, by = "sample_id")

# Create the plot for all samples collected on the specific date
ggplot(DateSpecific_merged, aes(x = sample_id, y = percent_value, fill = Family_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +
  theme_classic() +
  ylab("Relative Abundance") +
  xlab("Samples Collected on 2022-01-20") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


############################### Family Speed Pretest ###############################

# Filter metadata to keep only samples from the specific date
DateSpecific_metadata <- Metadata[Metadata$Date == "2022-02-17", ]

# Extract sample ID from row names, assuming rownames contain the sample IDs
DateSpecific_metadata$sample_id <- rownames(DateSpecific_metadata)

# Select only the necessary columns
DateSpecific_metadata <- DateSpecific_metadata[, c("sample_id", "Date")]

# Merge metadata with df_genus_combined
# Ensure that the ID column in df_family_combined matches the sample_id
DateSpecific_merged <- merge(DateSpecific_metadata, df_family_combined, by = "sample_id")

# Create the plot for all samples collected on the specific date
ggplot(DateSpecific_merged, aes(x = sample_id, y = percent_value, fill = Family_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +
  theme_classic() +
  ylab("Relative Abundance") +
  xlab("Samples Collected on 2021-12-07") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

################################### Family Filters and Pellets Pretest ################

# Filter metadata to keep only samples from the specific date
DateSpecific_metadata <- Metadata[Metadata$Date == "2022-04-21", ]

# Extract sample ID from row names, assuming rownames contain the sample IDs
DateSpecific_metadata$sample_id <- rownames(DateSpecific_metadata)

# Select only the necessary columns
DateSpecific_metadata <- DateSpecific_metadata[, c("sample_id", "Date")]

# Merge metadata with df_genus_combined
# Ensure that the ID column in df_family_combined matches the sample_id
DateSpecific_merged <- merge(DateSpecific_metadata, df_family_combined, by = "sample_id")

# Create the plot for all samples collected on the specific date
ggplot(DateSpecific_merged, aes(x = sample_id, y = percent_value, fill = Family_grouped)) +
  geom_col(position = "stack", color = "black", size = 0.5) +
  scale_fill_manual(values = my_palette) +
  theme_classic() +
  ylab("Relative Abundance") +
  xlab("Samples Collected on 2021-12-07") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))