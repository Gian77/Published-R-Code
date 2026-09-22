# ==============================================================================
# 01_import.R  -- build raw phyloseq objects for fungi (ITS) and prokaryotes (16S)
# Inputs : datasets/ITS/, datasets/16S/, mapping_ITS.txt, mapping_16s.txt
# Outputs: datasets/phyloseq/physeq_fungi_filt.rds, physeq_prokaryote_filt.rds
# ==============================================================================

# Locate the paper directory: this file is <ROOT>/code/<script>.R
ROOT <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(f)) dirname(dirname(normalizePath(f[1]))) else getwd()
})
source(file.path(ROOT, "functions", "00_setup.R"))

# ----------------------------------------------------------- A) ITS (fungi) ---
# OTU table (UPARSE, forward read R1)
otus_ITS <- read.delim(file.path(ITS_DIR, "otu_table_ITS_UPARSE_R1.txt"), row.names = 1)
otu_phy  <- otu_table(otus_ITS, taxa_are_rows = TRUE)

# metadata
meta_ITS <- read.delim(MAP_ITS, row.names = 1, header = TRUE, sep = "\t")
meta_phy <- sample_data(meta_ITS)

# taxonomy: CONSTAX2 consensus ranks (Kingdom..Species) + an extra "RDP" column
# (the RDP-classifier Kingdom call) used to select Fungi and flag the mock.
tax_cons <- read.delim(file.path(ITS_DIR, "constax_taxonomy.txt"),
                       header = TRUE, row.names = 1)
tax_rdp  <- read.delim(file.path(ITS_DIR, "rdp_taxonomy.txt"),
                       header = TRUE, row.names = 1)
tax_cons <- tax_cons[match(rownames(tax_rdp), rownames(tax_cons)), ]
tax_cons$RDP <- tax_rdp[[1 + 1]]          # RDP classifier Kingdom column ("Fungi"/"Mockk"/...)
tax_phy  <- tax_table(as.matrix(tax_cons))

seqs_ITS <- readDNAStringSet(file.path(ITS_DIR, "otus_R1.fasta"), format = "fasta",
                             seek.first.rec = TRUE, use.names = TRUE)

physeq_fungi <- phyloseq(otu_phy, meta_phy, tax_phy, seqs_ITS)
tax_table(physeq_fungi)[tax_table(physeq_fungi) == ""]   <- NA
tax_table(physeq_fungi)[is.na(tax_table(physeq_fungi))]  <- "Unclassified"

# Tag-switching was negligible; drop the mock sample (Ghcontrol1), filter rare
# OTUs (<10 reads total), and keep only RDP-confirmed Fungi.
otu_table(physeq_fungi) <- subset(otu_table(physeq_fungi), select = -c(Ghcontrol1))
otu_table(physeq_fungi) <- otu_table(physeq_fungi)[rowSums(otu_table(physeq_fungi)) >= 10, ]
physeq_fungi_filt <- subset_taxa(physeq_fungi, RDP == "Fungi")
message(sprintf("Fungi (filtered): %d OTUs x %d samples",
                ntaxa(physeq_fungi_filt), nsamples(physeq_fungi_filt)))

# --------------------------------------------------- B) 16S (prokaryotes) -----
otus_16s <- read.delim(file.path(BAC_DIR, "otu_table_16s_UPARSE.txt"), row.names = 1)
otu_phy  <- otu_table(otus_16s, taxa_are_rows = TRUE)

meta_16s <- read.delim(MAP_16S, row.names = 1, header = TRUE, sep = "\t")
meta_phy <- sample_data(meta_16s)

# taxonomy: RDP with a per-rank confidence cutoff of 0.7
tax_rdp16 <- read.delim(file.path(BAC_DIR, "taxonomy_16s_RDP.txt"), header = TRUE, row.names = 1)
keep_rank <- function(score, name) ifelse(score >= 0.7, as.character(name), NA)
taxonomy <- cbind(
  Kingdom = keep_rank(tax_rdp16$D_Score, tax_rdp16$Domain),
  Phylum  = keep_rank(tax_rdp16$P_Score, tax_rdp16$Phylum),
  Class   = keep_rank(tax_rdp16$C_Score, tax_rdp16$Class),
  Order   = keep_rank(tax_rdp16$O_Score, tax_rdp16$Order),
  Family  = keep_rank(tax_rdp16$F_Score, tax_rdp16$Family),
  Genus   = keep_rank(tax_rdp16$G_Score, tax_rdp16$Genus)
)
rownames(taxonomy) <- rownames(tax_rdp16)
tax_phy <- tax_table(as.matrix(taxonomy))

seqs_16s <- readDNAStringSet(file.path(BAC_DIR, "otus.fasta"), format = "fasta",
                             seek.first.rec = TRUE, use.names = TRUE)

physeq_prokaryote <- phyloseq(otu_phy, meta_phy, tax_phy, seqs_16s)
tax_table(physeq_prokaryote)[tax_table(physeq_prokaryote) == ""]  <- NA
tax_table(physeq_prokaryote)[is.na(tax_table(physeq_prokaryote))] <- "Unclassified"

# Drop unclassified-at-Kingdom and chloroplast/cyanobacteria, filter rare OTUs.
physeq_prokaryote <- subset_taxa(physeq_prokaryote, Kingdom != "Unclassified")
physeq_prokaryote <- subset_taxa(physeq_prokaryote, Phylum != "Cyanobacteria/Chloroplast")
otu_table(physeq_prokaryote) <- otu_table(physeq_prokaryote)[rowSums(otu_table(physeq_prokaryote)) >= 10, ]
physeq_prokaryote_filt <- physeq_prokaryote
message(sprintf("Prokaryotes (filtered): %d OTUs x %d samples",
                ntaxa(physeq_prokaryote_filt), nsamples(physeq_prokaryote_filt)))

# ------------------------------------------------------------------- save -----
save_physeq(physeq_fungi_filt,       "physeq_fungi_filt")
save_physeq(physeq_prokaryote_filt,  "physeq_prokaryote_filt")
