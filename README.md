# Tutorial how to get good accessions number and what can I do with them ? 

## 📋 Table of Contents

1. [Understand NCBI databases](#1-understand-ncbi-databases)
2. [Define your biological question](#2-define-your-biological-question)
3. [Master NCBI search syntax](#3-master-ncbi-search-syntax)
4. [Strategy A — Broad diversity sampling](#4-strategy-a--broad-diversity-sampling)
5. [Strategy B — Targeted strain retrieval](#5-strategy-b--targeted-strain-retrieval)
6. [Strategy C — Epidemiological outbreak sampling](#6-strategy-c--epidemiological-outbreak-sampling)
7. [Strategy D — Gene-level phylogeny (16S, ITS, MLST)](#7-strategy-d--gene-level-phylogeny-16s-its-mlst)
8. [Retrieve and download sequences](#8-retrieve-and-download-sequences)
9. [Add a root sequence (outgroup)](#9-add-a-root-sequence-outgroup)
10. [Align sequences with MAFFT](#10-align-sequences-with-mafft)
11. [Trim the alignment with ClipKIT](#11-trim-the-alignment-with-clipkit)
12. [Build the tree with RAxML](#12-build-the-tree-with-raxml)
13. [Visualize in R with ggtree](#13-visualize-in-r-with-ggtree)
14. [Quick reference — Organisms & recommended strategies](#14-quick-reference--organisms--recommended-strategies)

---
## 0. pre-require
Add you API-key
```Bash
API-KEY NCBI_API_KEY=855294aa6a5d0abb5c5286c6e93121a14908
```
set up micromamba environement 
```bash
micromamba create -n micro
micromamba activate micro
micromamba install -y -c conda-forge -c bioconda \
  entrez-direct \
  mafft \
  clipkit \
  raxml \
  r-base \
  bioconductor-ggtree \
  r-ggplot2 \
  r-ape
```
## 1. Understand NCBI databases

Before searching, you need to know **which NCBI database** to query.
Choosing the wrong one gives wrong or empty results.

| Database | Content | When to use it |
|---|---|---|
| `nucleotide` | DNA/RNA sequences | Genome, gene, 16S rRNA sequences |
| `protein` | Protein sequences | Enzyme, surface protein comparisons |
| `genome` | Complete assembled genomes | Whole-genome phylogenies |
| `assembly` | Genome assemblies + metadata | When you need strain-level metadata |
| `sra` | Raw sequencing reads | When you need raw Illumina/Nanopore data |

```bash
# Always specify the database with -db
esearch -db nucleotide -query "..."
esearch -db protein    -query "..."
esearch -db assembly   -query "..."
```

## 2. Define your biological question

**Before touching the terminal**, answer these three questions:

```
1. What organism am I studying?
   → Defines the [Organism] tag and the database to use

2. What is my phylogenetic goal?
   → Broad diversity? Outbreak tracing? Variant comparison? Gene evolution?

3. What sequence type do I need?
   → Complete genome? Single gene (16S, ITS, MLST)? Protein?
```

Your answers determine **which strategy to follow** (A, B, C or D below).

---

## 3. Master NCBI search syntax

NCBI uses **tagged fields** to filter results precisely.
Here are all the tags you will use in this pipeline:

### 3.1 — Organism and taxonomy filters

```bash
# Exact organism name (uses NCBI taxonomy)
"Salmonella enterica[Organism]"
"SARS-CoV-2[Organism]"
"Candida albicans[Organism]"
"Plasmodium falciparum[Organism]"

# Taxonomic group (all species within a genus)
"Salmonella[Organism]"
"Influenzavirus[Organism]"
"Aspergillus[Organism]"
```

### 3.2 — Sequence quality filters

```bash
# Only complete, fully assembled genomes
"complete genome[Title]"

# Exclude partial sequences
NOT "partial[Title]"

# Filter by sequence length (in base pairs)
"4000000:5500000[SLEN]"   # bacterial genomes ~4-5.5 Mb
"29000:32000[SLEN]"       # coronavirus genomes ~30 kb
"100:1600[SLEN]"          # 16S rRNA gene ~1500 bp
```

### 3.3 — Metadata filters

```bash
# Geographic origin
"France[Title]"
"Asia[Title]"

# Collection date range
"2020/01/01:2024/12/31[PDAT]"

# Host organism   /!\[Host] works only if it mensionned on the title ( sometime must be ignored or sited with [All feature] )
"Homo sapiens[Host]" 
"Gallus gallus[Host]"    # chicken
"Sus scrofa[Host]"       # pig

# Isolation source [Isolation sources] works only if it mensionned on the title ( sometime must be ignored or sited with **[All feature]** )
"food[All feature]"
"blood[All feature]"
"environment[All feature]"
```

### 3.4 — Combine filters with Boolean operators

```bash
# AND — both conditions must be true
"Salmonella enterica[Organism] AND complete genome[Title]"

# OR — at least one condition must be true
"Typhimurium[Title] OR Enteritidis[Title]"

# NOT — exclude a condition
"Salmonella enterica[Organism] NOT partial[Title]"

# Complex example
"Salmonella enterica[Organism] \
 AND complete genome[Title] \
 AND 2018/01/01:2024/01/01[PDAT] \
 AND Homo sapiens[Host] \
 NOT plasmid[Title]"
```

> ⚠️ **Always test your query first** by counting results before downloading:
> ```bash
> esearch -db nucleotide -query "YOUR QUERY HERE" | grep "Count" 
> ```

---

## 4. Strategy A — Broad diversity sampling

### "I want to see the global genetic diversity of my organism"

**Use case:** You have no prior hypothesis. You want an overview of all known
lineages, clades, or variants. Good for a first exploratory analysis.

**Recommended database:** `nucleotide` with `complete genome` filter

```bash
# Step 1 — Count available sequences
esearch -db nucleotide \
  -query "YOUR_ORGANISM[Organism] AND complete genome[Title]" \
  | grep "Count"

# Step 2 — Get a list of accession numbers to review
esearch -db nucleotide \
  -query "YOUR_ORGANISM[Organism] AND complete genome[Title]" \
  | efetch -format acc \
  > accession_list.txt

# Step 3 — Download up to 100 sequences
esearch -db nucleotide \
  -query "YOUR_ORGANISM[Organism] AND complete genome[Title]" \
  | efetch -format fasta -stop 100 \
  > MY_ORGANISM.fasta
```

**Organism-specific examples:**

```bash
# Salmonella enterica — complete chromosomes only, no plasmids
esearch -db nucleotide \
  -query "Salmonella enterica[Organism] \
          AND complete genome[Title] \
          NOT plasmid[Title]" \
  | efetch -format fasta -stop 100 \
  > salmonella.fasta

# SARS-CoV-2 — complete genomes, human host
esearch -db nucleotide \
  -query "SARS-CoV-2[Organism] \
          AND complete genome[Title] \
          AND Homo sapiens[Host]" \
  | efetch -format fasta -stop 100 \
  > sarscov2.fasta

# Candida albicans — fungal pathogen
esearch -db nucleotide \
  -query "Candida albicans[Organism] \
          AND complete genome[Title]" \
  | efetch -format fasta -stop 50 \
  > candida.fasta

# Mycobacterium tuberculosis — slow-growing bacterium
esearch -db nucleotide \
  -query "Mycobacterium tuberculosis[Organism] \
          AND complete genome[Title]" \
  | efetch -format fasta -stop 100 \
  > mtb.fasta
```

>  **tip:** NCBI returns sequences in arbitrary order by default.
> For true random sampling across diversity, use the `assembly` database
> and filter by `assembly_level` (see Strategy B).

---

## 5. Strategy B — Targeted strain retrieval

### "I have a list of specific strains I want to compare"

**Use case:** You already know which strains matter — reference strains,
outbreak isolates, or strains from a published paper.
You retrieve them **one by one by accession number**.

### 5.1 — Find accession numbers on the NCBI website

```
1. Go to https://www.ncbi.nlm.nih.gov/nucleotide/
2. Search: "Salmonella enterica Typhimurium LT2 complete genome"
3. Click the best matching result
4. The accession number is top-left: AE006468.2
5. Note the version number (the .2 after the dot — always use the latest)
```

### 5.2 — Find accession numbers programmatically

```bash
# Search and print accession numbers only
esearch -db nucleotide \
  -query "Salmonella enterica Typhimurium[Organism] complete genome[Title]" \
  | efetch -format acc

# Save to a file for review
esearch -db nucleotide \
  -query "Salmonella enterica Typhimurium[Organism] complete genome[Title]" \
  | efetch -format acc > typhimurium_accessions.txt

cat typhimurium_accessions.txt
```

### 5.3 — Retrieve a single strain by accession

```bash
# Generic command
esearch -db nucleotide -query "ACCESSION.VERSION" \
  | efetch -format fasta > ACCESSION.fasta

# Append to your main file
cat ACCESSION.fasta >> MY_ORGANISM.fasta
```

### 5.4 — Retrieve a batch of strains from a list

```bash
# If you have a file with one accession per line:
# accession_list.txt contains:
# AE006468.2
# AM933172.1
# AL513382.1
# ...

while read acc; do
  echo "Downloading $acc..."
  esearch -db nucleotide -query "$acc" \
    | efetch -format fasta >> MY_ORGANISM.fasta
  sleep 1  # be polite with NCBI servers
done < accession_list.txt
```

### 5.5 — Reference strains for common organisms

Here is a curated list of well-known reference strains with their accessions:

**Salmonella enterica**
| Strain | Serovar | Accession | Why use it |
|---|---|---|---|
| LT2 | Typhimurium | `AE006468.2` | Classic lab reference |
| P125109 | Enteritidis | `AM933172.1` | Major foodborne serovar |
| CT18 | Typhi | `AL513382.1` | Typhoid fever, good outgroup |
| SL1344 | Typhimurium | `FQ312003.1` | Virulence studies |
| Newport | Newport | `CP017251.1` | Multidrug-resistant |
| Heidelberg | Heidelberg | `CP001120.1` | Poultry outbreaks |
| Infantis | Infantis | `CP028537.1` | European poultry |

**SARS-CoV-2**
| Strain | Variant | Accession | Why use it |
|---|---|---|---|
| Wuhan-Hu-1 | Original | `NC_045512.2` | Universal reference |
| Alpha B.1.1.7 | Alpha | `MW761464.1` | UK variant |
| Delta B.1.617.2 | Delta | `MZ359841.1` | India variant |
| Omicron BA.1 | Omicron | `OL672836.1` | South Africa variant |

**Mycobacterium tuberculosis**
| Strain | Lineage | Accession | Why use it |
|---|---|---|---|
| H37Rv | Lineage 4 | `NC_000962.3` | Universal reference |
| CDC1551 | Lineage 4 | `AE000516.2` | Clinical isolate |
| BCG Pasteur | Vaccine | `AM408590.1` | Vaccine strain |

**Influenza A**
| Strain | Subtype | Accession | Why use it |
|---|---|---|---|
| A/Puerto Rico/8/1934 | H1N1 | `J02144.1` | Classic lab strain |
| A/California/07/2009 | H1N1pdm09 | `GQ117044.1` | 2009 pandemic |
| A/Hong Kong/1/1968 | H3N2 | `CY147325.1` | 1968 pandemic |

---

## 6. Strategy C — Epidemiological outbreak sampling

### "I want to trace a specific outbreak or transmission chain"

**Use case:** You have isolates from a known outbreak and want to place them
in a broader phylogenetic context. You combine:
- Your own outbreak sequences
- Background sequences from the same time period and region
- A reference strain as anchor

```bash
# Step 1 — Download background sequences from the same region and period
esearch -db nucleotide \
  -query "Salmonella enterica[Organism] \
          AND complete genome[Title] \
          AND France[Country] \
          AND 2022/01/01:2023/12/31[PDAT]" \
  | efetch -format fasta -stop 50 \
  > background.fasta

# Step 2 — Add your outbreak sequences (already in FASTA format)
cat my_outbreak_sequences.fasta >> background.fasta

# Step 3 — Add a reference strain as anchor
esearch -db nucleotide -query "AE006468.2" \
  | efetch -format fasta >> background.fasta

# Step 4 — Rename your file
mv background.fasta MY_ORGANISM.fasta
```

> **Important for outbreak analysis:**
> Make sure your outbreak sequences have **informative headers** in the FASTA file.
> A good header looks like: `>Sample_ID|Country|Date|Host`
> This metadata will appear on the tree leaves.

---

## 7. Strategy D — Gene-level phylogeny (16S, ITS, MLST)

### "I want to compare a specific gene, not the whole genome"

**Use case:** Whole genomes are too large, too few, or unavailable.
You use a **marker gene** instead:
- **16S rRNA** → universal bacterial phylogeny
- **ITS** → fungal phylogeny
- **rbcL / matK** → plant phylogeny
- **COI** → animal/metazoan phylogeny
- **MLST genes** → fine-scale bacterial epidemiology

### 7.1 — 16S rRNA (bacteria)

```bash
# Search 16S sequences — note the length filter (1400-1600 bp)
esearch -db nucleotide \
  -query "Salmonella enterica[Organism] \
          AND 16S ribosomal RNA[Title] \
          AND 1400:1600[SLEN]" \
  | efetch -format fasta -stop 100 \
  > salmonella_16S.fasta
```

### 7.2 — ITS region (fungi)

```bash
# ITS1-5.8S-ITS2 region for fungal identification
esearch -db nucleotide \
  -query "Candida albicans[Organism] \
          AND internal transcribed spacer[Title] \
          AND 400:800[SLEN]" \
  | efetch -format fasta -stop 100 \
  > candida_ITS.fasta
```

### 7.3 — COI gene (animals / metabarcoding)

```bash
# Cytochrome oxidase I — the universal animal barcode
esearch -db nucleotide \
  -query "Aedes aegypti[Organism] \
          AND cytochrome oxidase subunit I[Title] \
          AND 600:700[SLEN]" \
  | efetch -format fasta -stop 100 \
  > aedes_COI.fasta
```

### 7.4 — MLST genes (fine-scale bacterial epidemiology)

```bash
# Example: housekeeping gene aroC for Salmonella MLST
esearch -db nucleotide \
  -query "Salmonella enterica[Organism] \
          AND aroC[Gene] \
          AND 400:600[SLEN]" \
  | efetch -format fasta -stop 100 \
  > salmonella_aroC.fasta
```

> 💡 **Which strategy for which question?**
>
> | Question | Strategy | Sequence type |
> |---|---|---|
> | Global diversity overview | A | Complete genome |
> | Compare known strains | B | Complete genome |
> | Trace an outbreak | C | Complete genome |
> | Universal bacterial ID | D | 16S rRNA |
> | Universal fungal ID | D | ITS |
> | Universal animal ID | D | COI |
> | Fine-scale epidemiology | D | MLST genes |

---

## 8. Retrieve and download sequences

Once your strategy is defined, here is the **complete download workflow**:

```bash
# Step 1 — Always count first
esearch -db nucleotide -query "YOUR COMPLETE QUERY" | grep "Count"

# Step 2 — Preview accession numbers (no download yet)
esearch -db nucleotide -query "YOUR COMPLETE QUERY" \
  | efetch -format acc | head -20

# Step 3 — Save the full accession list
esearch -db nucleotide -query "YOUR COMPLETE QUERY" \
  | efetch -format acc > accession_list.txt

# Step 4 — Download sequences in FASTA format
esearch -db nucleotide -query "YOUR COMPLETE QUERY" \
  | efetch -format fasta -stop 100 > MY_ORGANISM.fasta

# Step 5 — Verify the download
grep ">" MY_ORGANISM.fasta | wc -l       # count sequences
grep ">" MY_ORGANISM.fasta | head -10    # preview headers
```

> ⚠️ **NCBI rate limits:**
> - Without API key: max 3 requests/second
> - With API key: max 10 requests/second
>
> ```bash
> # Set your API key (get it free at ncbi.nlm.nih.gov/account)
> export NCBI_API_KEY=your_key_here
> ```

---

## 9. Add a root sequence (outgroup)

An **outgroup** is a distantly related sequence used to **root the tree**.
Without it, the tree has no direction — you cannot tell which lineages
are ancestral and which are derived.

```bash
# Rule: the outgroup should be related but clearly distinct
# from all your ingroup sequences

# For Salmonella enterica → use Escherichia coli K-12
esearch -db nucleotide -query "U00096.3" \
  | efetch -format fasta > outgroup.fasta
cat outgroup.fasta >> MY_ORGANISM.fasta

# For SARS-CoV-2 → use SARS-CoV-1 (2003 epidemic)
esearch -db nucleotide -query "AY274119.3" \
  | efetch -format fasta > outgroup.fasta
cat outgroup.fasta >> MY_ORGANISM.fasta

# For Candida albicans → use Saccharomyces cerevisiae
esearch -db nucleotide -query "NC_001133.9" \
  | efetch -format fasta > outgroup.fasta
cat outgroup.fasta >> MY_ORGANISM.fasta

# For Mycobacterium tuberculosis → use Mycobacterium smegmatis
esearch -db nucleotide -query "CP000480.1" \
  | efetch -format fasta > outgroup.fasta
cat outgroup.fasta >> MY_ORGANISM.fasta
```

> 💡 **How to choose an outgroup:**
> - It must be in the **same family or order** as your organism
> - It must be **clearly outside** your ingroup
> - It should have a **complete, high-quality genome**
> - When in doubt, check published phylogenies for your organism

---

## 10. Align sequences with MAFFT

Raw sequences are unaligned — each sequence starts at a different position.
MAFFT lines them up **column by column** so homologous positions are compared.

```bash
# Standard alignment — good for most cases
mafft --auto --thread 8 MY_ORGANISM.fasta > MY_ORGANISM.mafft.fasta

# For large whole genomes (>100 sequences, >1 Mb each) — faster
mafft --retree 1 --thread 16 MY_ORGANISM.fasta > MY_ORGANISM.mafft.fasta

# For highly divergent sequences — more accurate
mafft --localpair --maxiterate 1000 --thread 8 MY_ORGANISM.fasta \
  > MY_ORGANISM.mafft.fasta
```

> ⏱️ **Expected time:**
> - 100 bacterial genomes (~4 Mb each): 15–45 minutes
> - 100 viral genomes (~30 kb each): 1–3 minutes
> - 100 16S sequences (~1.5 kb each): < 1 minute

---

## 11. Trim the alignment with ClipKIT

Poorly aligned regions and gaps add noise to the tree.
ClipKIT removes them automatically.

```bash
# Default trimming — removes gappy columns
clipkit MY_ORGANISM.mafft.fasta -o MY_ORGANISM.trimmed.fasta

# Check how many positions were kept
grep -v ">" MY_ORGANISM.mafft.fasta | head -1 | wc -c    # before
grep -v ">" MY_ORGANISM.trimmed.fasta | head -1 | wc -c  # after
```

---

## 12. Build the tree with RAxML

RAxML uses **Maximum Likelihood** to find the most probable evolutionary tree.

```bash
# Standard run with 100 bootstrap replicates
raxmlHPC -f a \
  -m GTRGAMMA \
  -p 12345 \
  -x 12345 \
  -# 100 \
  -s MY_ORGANISM.trimmed.fasta \
  -n MY_ORGANISM_TREE \
  -T 8

# The final tree is in:
# RAxML_bipartitions.MY_ORGANISM_TREE
```

> 💡 **Model selection:**
> - `GTRGAMMA` → standard for DNA sequences (bacteria, viruses)
> - `PROTGAMMAAUTO` → for protein sequences
> - `BINCAT` → for SNP-only alignments

---

## 13. Visualize in R with ggtree

```r
# Install packages (once)
install.packages("BiocManager")
BiocManager::install("ggtree")
install.packages("ggplot2")

library(ggtree)
library(ggplot2)

# Load the tree
tree <- read.tree("RAxML_bipartitions.MY_ORGANISM_TREE")

# Root the tree on the outgroup
# Replace "outgroup_sequence_name" with the actual FASTA header of your outgroup
tree <- root(tree, outgroup = "outgroup_sequence_name", resolve.root = TRUE)

# Basic tree plot
ggtree(tree) +
  geom_tiplab(size = 2) +
  geom_nodelab(aes(label = label), size = 2, color = "red") +
  theme_tree2() +
  ggtitle("Phylogenetic tree — MY ORGANISM")

# Save the figure
ggsave("phylogenetic_tree.pdf", width = 12, height = 16)
```

---

## 14. Quick reference — Organisms & recommended strategies

| Organism | Strategy | Database | Sequence type | Outgroup |
|---|---|---|---|---|
| Salmonella enterica | A or B | nucleotide | Complete genome | E. coli K-12 `U00096.3` |
| SARS-CoV-2 | A or C | nucleotide | Complete genome | SARS-CoV-1 `AY274119.3` |
| Mycobacterium tuberculosis | B or C | nucleotide | Complete genome | M. smegmatis `CP000480.1` |
| Candida albicans | A or D | nucleotide | Complete genome or ITS | S. cerevisiae `NC_001133.9` |
| Influenza A | A or B | nucleotide | Complete genome or HA gene | Influenza B `CY115151.1` |
| Plasmodium falciparum | B | nucleotide | Complete genome | P. vivax `CM000450.1` |
| E. coli (pathogenic) | A or C | nucleotide | Complete genome | Salmonella LT2 `AE006468.2` |
| Any bacterium (ID only) | D | nucleotide | 16S rRNA | Archaea 16S |
| Any fungus (ID only) | D | nucleotide | ITS | Distant fungal genus |
| Any animal (barcoding) | D | nucleotide | COI | Distant animal phylum |

---

## ✅ Full pipeline summary

```bash
# 1. Search and count
esearch -db nucleotide -query "YOUR_ORGANISM[Organism] AND complete genome[Title]" \
  | grep "Count"

# 2. Get accession list
esearch -db nucleotide -query "YOUR_ORGANISM[Organism] AND complete genome[Title]" \
  | efetch -format acc > accession_list.txt

# 3. Download sequences
esearch -db nucleotide -query "YOUR_ORGANISM[Organism] AND complete genome[Title]" \
  | efetch -format fasta -stop 100 > MY_ORGANISM.fasta

# 4. Add reference strains (repeat for each)
esearch -db nucleotide -query "ACCESSION.VERSION" \
  | efetch -format fasta >> MY_ORGANISM.fasta

# 5. Add outgroup
esearch -db nucleotide -query "OUTGROUP_ACCESSION" \
  | efetch -format fasta >> MY_ORGANISM.fasta

# 6. Align
mafft --auto --thread 8 MY_ORGANISM.fasta > MY_ORGANISM.mafft.fasta

# 7. Trim
clipkit MY_ORGANISM.mafft.fasta -o MY_ORGANISM.trimmed.fasta

# 8. Build tree
raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -# 100 \
  -s MY_ORGANISM.trimmed.fasta -n MY_ORGANISM_TREE -T 8

# 9. Visualize → open R and run the ggtree script
```
