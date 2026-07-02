# 📘 Code for Kelley et al. Arabidopsis Drought SynCom

## **Title:** "Investigating the Influence of Drought on the Assembly of a Designed Root Microbial Community in Arabidopsis"

## 📖 Citation (unpublished)

**Authors:** Brittni R. Kelley 1; Gian Maria Niccolò **Benucci** 2,3; Marlana A. DeClaire 3,4;
Brad Day 5; Sarah Lebeis 1,2,3,4

1 Plant Resilience Institute, Michigan State University<br>
2 Great Lakes Bioenergy Research Center, Michigan State University<br>
3 Department of Plant, Soil, and Microbial Sciences, Michigan State University, East Lansing, MI 48823, USA<br>
4 Microbiology and Molecular Genetics, Michigan State University<br>
5 Microbiology, University of Tennessee-Knoxville

Corresponding author: benucci@msu.edu

## **Journal:** to be selected

**Status:** In preparation.

---

## 📂 Repository Structure

`R/` holds the numbered analysis pipeline (see `CLAUDE.md` for script order); `code/` holds
legacy monolithic scripts kept for reference only.

```
├── code
│   ├── Rcode_SynComDrought26.R
│   └── Rcode_SynComDrought.R
├── datasets
│   ├── asv_233bp.fasta
│   ├── asv_233bp.sintax
│   ├── asv_233bp_to97otus.fasta
│   ├── asv_233bp_to97otus.sintax
│   ├── asvtable_UNOISE_233bp.txt
│   ├── file_mapping.txt
│   ├── metadata.csv
│   ├── otus_233bp.fasta
│   ├── otus_233bp.sintax
│   ├── otutable_233bp_closedRef.txt
│   ├── otutable_asv_233bp_to97otus.txt
│   ├── otutable_SWARM_233bp.txt
│   ├── otutable_UPARSE_233bp.txt
│   ├── swarms_233bp.fasta
│   └── swarms_233bp.sintax
├── functions
├── misc
├── R
│   ├── 00_setup_evnviron.R
│   └── 01_import_data.R
├── renv
├── results
│   ├── bpLength-faceted.pdf
│   ├── bpLength-log2.pdf
│   └── bpLength.pdf
├── skills
├── CLAUDE.md
├── LICENSE
└── README.md
```

## 🚀 Clone and Use
```
# clone the main repo
git clone https://github.com/Gian77/Published-R-Code.git

# get only the desired subrepo
svn export https://github.com/Gian77/Published-R-Code/tree/master/Kelley_etal_2024_ArabidopsisDrought

cd Kelley_etal_2024_ArabidopsisDrought
```
