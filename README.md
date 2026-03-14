# 🧬 Phylogenetic Analysis Pipeline — Any Organism

## 📋 Table of Contents

1. [Prerequisites](#1-prerequisites)
2. [Setup your environment](#2-setup-your-environment)
3. [Search & Download sequences from NCBI](#3-search--download-sequences-from-ncbi)
4. [Add specific reference strains](#4-add-specific-reference-strains)
5. [Align sequences with MAFFT](#5-align-sequences-with-mafft)
6. [Trim the alignment with ClipKIT](#6-trim-the-alignment-with-clipkit)
7. [Build the phylogenetic tree with RAxML](#7-build-the-phylogenetic-tree-with-raxml)
8. [Visualize the tree in R](#8-visualize-the-tree-in-r)
9. [Adapt this pipeline to your organism](#9-adapt-this-pipeline-to-your-organism)

---

## 1. Prerequisites

Make sure the following tools are installed.
We use **micromamba** (a lightweight conda alternative) to manage them.

```bash
# Install all tools at once
micromamba install mafft clipkit raxml
