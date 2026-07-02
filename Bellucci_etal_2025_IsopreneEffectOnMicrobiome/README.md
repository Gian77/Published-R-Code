# 📘 Code for Bellucci et al. Isoprene-Emitting Tobacco Microbiome

## **Title:** "Isoprene-emitting Transgenic Tobacco Shapes Root Microbiome and Enhances Growth of Co-cultivated Non-emitting Plants."

## 📖 Citation (unpublished)

**Authors:** Bellucci M, Mostofa M, **Benucci GMN**, Kabir A, Khan I, Lombardi M, Locato V,
Bonito G, Loreto F, Sharkey T

Corresponding author: benucci@msu.edu

## **Journal:** to be selected

**Status:** Submitted.

_NOTE_ The ITS and 16S datasets are cleaned of contaminant OTUs (using the `decontam` R
package), non-target taxa (e.g. mitochondria), and control samples. Taxonomy and metadata have
been cleaned for ease of use as well.

---

## 📂 Repository Structure

```
├── code
│   └── Rcode_IsoprenePlantMicrobiome.R
├── datasets
│   ├── phyloseq_16S.rds
│   └── phyloseq_ITS.rds
├── functions
├── misc
├── LICENSE
└── README.md
```

## 🚀 Clone and Use
```
# clone the main repo
git clone https://github.com/Gian77/Published-R-Code.git

# get only the desired subrepo
svn export https://github.com/Gian77/Published-R-Code/tree/master/Bellucci_etal_2025_IsopreneEffectOnMicrobiome

cd Bellucci_etal_2025_IsopreneEffectOnMicrobiome
```
