## Code used in:

Isoprene-emitting Transgenic Tobacco Shapes Root Microbiome and Enhances Growth of Co-cultivated Non-emitting Plants. Bellucci M, Mostofa M, Benucci GMN, Kabir A, Khan I, Lombardi M, Locato V, Bonito G, Loreto F, Sharkey T (submitted)

To reprodu e the analysis in R. Download the `.rds` and load it as:

```
# Load the RDS file back into an R object
phyloseq_ITS <- readRDS("phyloseq_ITS.rds")
phyloseq_16S <- readRDS("phyloseq_16S.rds")
```
_NOTE_ The ITS and 16S dataset are clean out of contaminant OTUs (using `decontam` R package), non-target taxa (e.g. mitochondria), and control samples. The taxonomy and the metadata have been cleaned for ease to use as well. 


The follow the provided `.R` code.
--

