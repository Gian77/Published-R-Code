# 📘 Code for Benucci et al. (2019) Front. Microbiol.

## 📖 Citation

**Benucci GMN**, Burnard D, Shepherd LD, Bonito G and Munkacsi AB. Evidence for 
Co-evolutionary History of Early Diverging Lycopodiaceae Plants With Fungi. Front 
Microbiol. 2020 Jan 15;10:2944. doi: 10.3389/fmicb.2019.02944. PMID: 32010072; 
PMCID: PMC6974469.

Corresponding author: Andrew B. Munkacsi andrew.munkacsi@vuw.ac.nz

DOI: [10.3389/fmicb.2019.02944](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2019.02944/full) <br>
URL: https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2019.02944/full

---

## 📂 Repository Structure

```
├── code
│   └── Rcode_LycopodsNZ.R
├── datasets
│   ├── mapping_16s.txt
│   ├── mapping_18S.txt
│   ├── mapping_ITS.txt
│   ├── mapping_SSU.txt
│   ├── otus_16S.fasta
│   ├── otus_R1_18S.fasta
│   ├── otus_R1_ITS.fasta
│   ├── otus_R1_SSU.fasta
│   ├── otu_table_16S_UPARSE.txt
│   ├── otu_table_18S_UPARSE_R1.txt
│   ├── otu_table_ITS_UPARSE_R1.txt
│   ├── otu_table_SSU_UPARSE_R1.txt
│   ├── taxonomy_16S_RDP.txt
│   ├── taxonomy_16S_SILVA.txt
│   ├── taxonomy_18S_SILVA.csv
│   ├── taxonomy_ITS_consensus.txt
│   ├── taxonomy_ITS_RDP.txt
│   └── taxonomy_SSU_SILVA.csv
├── functions
└── README.md
```

## 🚀 Clone and Use
```
# clone the main repo
git clone https://github.com/Gian77/Published-R-Code.git

# get only the desired subrepo
svn export https://github.com/Gian77/Published-R-Code/tree/master/Benucci_etal_2019_Lycopods

cd Benucci_etal_2019_Lycopods
```