# Phylogenetic Analysis of Kinesin Motor and Kinesin-14 Family Proteins

A complete end-to-end pipeline for phylogenetic analysis of two kinesin protein families — from remote sequence retrieval and BLAST search through multiple sequence alignment, model selection, maximum-likelihood tree inference, bootstrap support estimation, and ancestral sequence reconstruction.

---

## Table of Contents

- [Background](#background)
- [Pipeline Overview](#pipeline-overview)
- [Dependencies and Requirements](#dependencies-and-requirements)
- [Installation](#installation)
- [Project Structure](#project-structure)
- [Methods](#methods)
- [Step-by-Step Usage](#step-by-step-usage)
- [Output Files](#output-files)
- [Expected Output Details](#expected-output-details)
- [Key Biological Interpretation](#key-biological-interpretation)
- [References](#references)

---

## Background

This project performs a standard computational phylogenetics workflow applied to two protein families from the Duke Kinesin Resource:

- **Kinesin motor proteins** (Kinesin-1 family): conventional N-terminal motors, plus-end directed, involved in vesicle and organelle transport. The canonical query sequence is human **KIF5B** (NCBI: `NP_004512`).
- **Kinesin-14 family proteins**: C-terminal motors, minus-end directed, involved in spindle assembly and chromosome segregation. Key members include *Drosophila* Ncd, *S. cerevisiae* Kar3, *S. pombe* Pkl1, *Xenopus* XCTK2, and human KIFC1/KIFC2/KIFC3. The canonical query sequence is human **KIFC1** (NCBI: `NP_002254`), also known as HSET.

Two separate BLASTp searches are performed — one per query — to retrieve diverse cross-species homologs. The hits are merged and processed through a unified pipeline to produce a tree showing both families as distinct, well-supported clades.

---

## Pipeline Overview

```
Query A: KIF5B (Kinesin motor, NP_004512)       Query B: KIFC1 (Kinesin-14, NP_002254)
                    │                                              │
                    ▼                                              ▼
         BLASTp vs PDB (top 15 hits)                  BLASTp vs PDB (top 15 hits)
                    │                                              │
                    └──────────────────┬───────────────────────────┘
                                       ▼
                          Merge + Deduplicate Sequences
                                       │
                                       ▼
                      Multiple Sequence Alignment (MAFFT)
                                       │
                                       ▼
                         Trim Alignment (gappy ends)
                                       │
                                       ▼
                        Remove Duplicate Sequences
                                       │
                                       ▼
                       Model Selection (ModelTest-NG)
                                       │
                                       ▼
                     ML Tree Inference (RAxML-NG, LG+I+G4)
                                       │
                                       ▼
                      Bootstrap Analysis (100 replicates)
                                       │
                                       ▼
                     Support Mapping → Final Annotated Tree
                                       │
                                       ▼
                     Ancestral Sequence Reconstruction
```

---

## Dependencies and Requirements

| Tool / Library | Version Tested | Purpose                                      |
|----------------|----------------|----------------------------------------------|
| Python         | >= 3.8         | Core language                                |
| BioPython      | >= 1.79        | Sequence I/O, BLAST, alignment, tree parsing |
| MAFFT          | >= 7.0         | Multiple sequence alignment                  |
| ModelTest-NG   | >= 0.1.7       | Substitution model selection                 |
| RAxML-NG       | >= 1.2.0       | Maximum-likelihood tree inference            |
| matplotlib     | >= 3.0         | Tree visualization                           |

> **Note:** This pipeline was developed and tested on Google Colab. RAxML-NG and ModelTest-NG are installed via direct binary download — no conda required.

---

## Installation

### Google Colab 

```python
# Install BioPython and MAFFT
!pip install biopython -q
!apt-get install -qq -y mafft

# Install RAxML-NG and ModelTest-NG via direct binary download
!wget -q <raxml-ng-release-url> && unzip && cp to /usr/local/bin
!wget -q <modeltest-ng-release-url> && tar -xzf && cp to /usr/local/bin
```

> **Important:** Do NOT use `condacolab` — it fails on current Colab environments due to a SHA256 checksum mismatch with the Miniconda installer.

### Local Installation

```bash
conda create -n phylo python=3.9
conda activate phylo
conda install -c conda-forge -c bioconda biopython mafft modeltest-ng raxml-ng matplotlib
```

---

## Project Structure

```
phylogenetic-analysis/
│
├── data/
│   ├── KIF5B_Hs_Kinesin1.fasta          # Query A — Kinesin motor (KIF5B)
│   ├── KIFC1_Hs_Kinesin14.fasta         # Query B — Kinesin-14 (KIFC1/HSET)
│   ├── sequences.fasta                  # All BLAST hits merged with queries
│   ├── aligned.fasta                    # MAFFT alignment output
│   ├── aligned_trimmed.fasta            # Trimmed alignment (gappy ends removed)
│   ├── clear_aligned_trimmed.fasta      # Deduplicated final alignment
│   ├── tree_ML.png                      # Initial ML tree visualization
│   ├── tree_bootstrap.png               # Bootstrap-annotated tree visualization
│   └── tree_ancestral.png               # Ancestral reconstruction tree visualization
│
├── T1.raxml.bestTree                    # Best ML tree (Newick format)
├── T1.raxml.rba                         # Binary alignment for bootstrapping
├── T1.raxml.mlTrees                     # All ML trees from multiple starts
├── T2.raxml.bootstraps                  # 100 bootstrap replicate trees
├── T3.raxml.support                     # ML tree annotated with bootstrap support
├── ASR.raxml.ancestralTree              # Ancestral reconstruction tree
├── ASR.raxml.ancestralProbs             # Per-node ancestral amino acid probabilities
│
└── README.md
```

---

## Methods

Two query sequences (KIF5B and KIFC1) were fetched from NCBI, and two BLASTp searches against the PDB database retrieved up to 15 cross-species hits per query. Hits were merged, deduplicated, and aligned with **MAFFT**. The alignment was end-trimmed and cleaned of duplicates. **ModelTest-NG** selected **LG+I+G4** as the best substitution model by BIC. A maximum-likelihood tree was built with **RAxML-NG**, branch support assessed via 100 bootstrap replicates, and marginal ancestral amino acid states reconstructed at all internal nodes.

---

## Step-by-Step Usage

Open the notebook and run cells sequentially. Each step is self-contained and annotated.

---

### Step 0A — Install BioPython and MAFFT

Installs BioPython for sequence handling and MAFFT for multiple sequence alignment.

```python
# Install BioPython and MAFFT, then verify MAFFT version
!pip install biopython -q
!apt-get install -qq -y mafft
```

---

### Step 0B — Install ModelTest-NG and RAxML-NG

Downloads pre-compiled binaries for ModelTest-NG and RAxML-NG directly from GitHub releases and places them on the system path.

```python
# Download RAxML-NG and ModelTest-NG binaries and copy to /usr/local/bin
!wget raxml-ng → unzip → cp /usr/local/bin
!wget modeltest-ng → tar -xzf → cp /usr/local/bin
```

---

### Step 0C — Verify Installation

Confirms both tools are on the system path before proceeding.

```python
# Check both tools are accessible and print their versions
!which modeltest-ng && modeltest-ng --version
!which raxml-ng    && raxml-ng --version
```

---

### Step 1 — Retrieve Both Query Sequences

Fetches KIF5B (`NP_004512`) and KIFC1 (`NP_002254`) from NCBI Protein using `Entrez.efetch` and saves each as a FASTA file.

```python
# Fetch KIF5B (Kinesin-1) and KIFC1 (Kinesin-14) from NCBI
# Save each to data/ as a named FASTA file
Entrez.efetch(db="protein", id="NP_004512") → save KIF5B_Hs_Kinesin1.fasta
Entrez.efetch(db="protein", id="NP_002254") → save KIFC1_Hs_Kinesin14.fasta
```

---

### Step 2 — Run Two Separate BLAST Searches

Performs two BLASTp searches against PDB with no organism filter — one per query — to capture cross-species kinesin motor and kinesin-14 homologs. Each search may take 3–5 minutes.

```python
# BLASTp search 1: KIF5B vs PDB (top 30 hits) → save XML
NCBIWWW.qblast("blastp", "pdb", KIF5B.seq, hitlist_size=30) → blast_kinesin_motor.xml

# BLASTp search 2: KIFC1 vs PDB (top 30 hits) → save XML
NCBIWWW.qblast("blastp", "pdb", KIFC1.seq, hitlist_size=30) → blast_kinesin14.xml
```

---

### Step 3 — Parse Both BLAST Results

Reads the saved XML files into structured `blast_record` objects using `NCBIXML`.

```python
# Parse both XML result files into blast_record objects
NCBIXML.read("blast_kinesin_motor.xml") → blast_kinesin
NCBIXML.read("blast_kinesin14.xml")     → blast_k14
```

---

### Step 4 — View BLAST Hits for Both Searches

Iterates over all alignments from both searches and prints Hit ID, alignment length, E-value, and percent identity.

```python
# Print Hit ID, alignment length, E-value, PID% for each hit
for aln in blast_record.alignments:
    print(hit_id, align_length, e_value, percent_identity)
```

---

### Step 5 — Filter Hits by Percent Identity

Applies a percent identity cutoff (`PIDcut = 1.00`) and reports sequence length, alignment length, E-value, PID%, and query coverage for each passing hit.

```python
# Keep hits where PID <= PIDcut (1.00 = accept all)
# Report accession, seq length, alignment length, e-value, PID%, coverage%
if (identities / align_length) <= PIDcut → print and collect hit
```

---

### Step 6 — Download Top Hits and Merge into One FASTA

Downloads the top 15 hits from each search, adds both query sequences as anchors, deduplicates by exact sequence, and writes the combined set to `data/sequences.fasta`.

```python
# Fetch top 15 hits per search via Entrez.efetch
# Prepend both query sequences, deduplicate by sequence string
# Write all unique sequences → data/sequences.fasta
```

---

### Step 7 — Multiple Sequence Alignment (MAFFT)

Aligns all sequences using MAFFT auto strategy. Output saved to `data/aligned.fasta`.

```python
# Align all sequences with MAFFT auto strategy
!mafft --auto data/sequences.fasta > data/aligned.fasta
```

---

### Step 8 — Trim Alignment

Removes terminal gap-containing columns from both ends by scanning left-to-right for the first ungapped column and right-to-left for the last, then slicing to that range.

```python
# Find first ungapped column from left  → position1
# Find last  ungapped column from right → position2
# Slice alignment[:, position1:position2] → aligned_trimmed.fasta
```

---

### Step 9 — Compute Pairwise Distance Matrix

Computes a pairwise identity-based distance matrix from the trimmed alignment to give a numerical overview of sequence divergence before tree building.

```python
# Compute pairwise distances using identity model
calculator = DistanceCalculator("identity")
dm = calculator.get_distance(trimmed_alignment)
```

---

### Step 10 — Remove Duplicate Sequences

Runs `sequence_cleaner()` to remove exact duplicate sequences, merging their IDs with an underscore. Output saved as `data/clear_aligned_trimmed.fasta`.

```python
# Remove exact duplicate sequences, merge duplicate IDs with underscore
# Input: aligned_trimmed.fasta → Output: clear_aligned_trimmed.fasta
sequence_cleaner("data/aligned_trimmed.fasta")
```

---

### Step 11 — Model Selection (ModelTest-NG)

Runs ModelTest-NG on the cleaned alignment to find the best-fit amino acid substitution model by BIC. Expected result: **LG+I+G4**.

```python
# Run ModelTest-NG, select best model by BIC
!modeltest-ng -i clear_aligned_trimmed.fasta -d aa
MODEL = "LG+I+G4"   # update based on BIC output above
```

---

### Step 12 — Build Maximum-Likelihood Tree (RAxML-NG)

Infers the best-scoring ML tree under the selected model. Produces `T1.raxml.bestTree` (Newick tree) and `T1.raxml.rba` (binary alignment for bootstrapping).

```python
# Infer best ML tree under selected model, prefix T1
!raxml-ng --msa clear_aligned_trimmed.fasta --model LG+I+G4 --prefix T1 --seed 12345
```

---

### Step 13 — Visualize Initial Tree

Reads `T1.raxml.bestTree` and renders it with matplotlib. KIF5B highlighted in blue, KIFC1 in red.

```python
# Read T1.raxml.bestTree, draw with Phylo, highlight KIF5B (blue) and KIFC1 (red)
tree = Phylo.read("T1.raxml.bestTree", "newick")
Phylo.draw(tree) → save data/tree_ML.png
```

---

### Step 14 — Bootstrap Analysis

Runs 100 bootstrap replicates using `T1.raxml.rba`. Bootstrap trees saved to `T2.raxml.bootstraps`.

```python
# Run 100 bootstrap replicates on binary alignment from Step 12
!raxml-ng --bootstrap --msa T1.raxml.rba --model LG+I+G4 --prefix T2 --bs-trees 100
```

---

### Step 15 — Map Bootstrap Support Values

Maps the 100 bootstrap trees onto the best ML tree. Annotated tree saved to `T3.raxml.support`.

```python
# Map bootstrap support from T2 onto best ML tree from T1 → T3.raxml.support
!raxml-ng --support --tree T1.raxml.bestTree --bs-trees T2.raxml.bootstraps --prefix T3
```

---

### Step 16 — Visualize Bootstrap Support Tree

Reads `T3.raxml.support` and renders the tree with bootstrap percentages at internal nodes.

```python
# Read T3.raxml.support, draw tree with node support values
# KIF5B (blue), KIFC1 (red) → save data/tree_bootstrap.png
tree = Phylo.read("T3.raxml.support", "newick")
Phylo.draw(tree) → save data/tree_bootstrap.png
```

---

### Step 17 — Ancestral Sequence Reconstruction

Runs RAxML-NG in `--ancestral` mode to reconstruct marginal ancestral amino acid states at every internal node.

```python
# Reconstruct ancestral amino acid states at all internal nodes
!raxml-ng --ancestral --msa clear_aligned_trimmed.fasta --tree T1.raxml.bestTree --model LG+I+G4 --prefix ASR
```

---

### Step 18 — Visualize Ancestral Tree

Reads `ASR.raxml.ancestralTree` and renders it with all internal nodes labelled for cross-referencing with the ancestral probability file.

```python
# Read ASR.raxml.ancestralTree, draw with labelled internal nodes
# Internal node labels reference ASR.raxml.ancestralProbs
tree = Phylo.read("ASR.raxml.ancestralTree", "newick")
Phylo.draw(tree) → save data/tree_ancestral.png
```

---

### Download All Results

```python
# Download all output files to local machine
files.download("data/sequences.fasta")
files.download("T1.raxml.bestTree")
files.download("T3.raxml.support")
files.download("ASR.raxml.ancestralProbs")
# ... and all other output files
```

---

## Output Files

| File | Description |
|------|-------------|
| `data/KIF5B_Hs_Kinesin1.fasta` | Query A — Human KIF5B, Kinesin-1 motor (963 aa) |
| `data/KIFC1_Hs_Kinesin14.fasta` | Query B — Human KIFC1/HSET, Kinesin-14 (673 aa) |
| `data/sequences.fasta` | All BLAST hits merged with both query sequences |
| `data/aligned.fasta` | MAFFT multiple sequence alignment (raw) |
| `data/aligned_trimmed.fasta` | Alignment after removing gappy terminal columns |
| `data/clear_aligned_trimmed.fasta` | Final deduplicated alignment used for tree building |
| `T1.raxml.bestTree` | Best-scoring maximum-likelihood tree (Newick format) |
| `T1.raxml.rba` | Binary compressed alignment used for bootstrapping |
| `T1.raxml.mlTrees` | All ML trees inferred from multiple starting points |
| `T2.raxml.bootstraps` | 100 bootstrap replicate trees |
| `T3.raxml.support` | Best ML tree annotated with bootstrap support values |
| `ASR.raxml.ancestralTree` | Tree with labelled internal nodes for ancestral states |
| `ASR.raxml.ancestralProbs` | Per-node marginal amino acid probability matrix |
| `data/tree_ML.png` | Visualization of best ML tree (no bootstrap) |
| `data/tree_bootstrap.png` | Visualization of ML tree with bootstrap support |
| `data/tree_ancestral.png` | Visualization of ancestral reconstruction tree |

---

## Expected Output Details

### Step 1 — Sequence Retrieval
```
[KIF5B_Hs_Kinesin1]   length=963   saved → data/KIF5B_Hs_Kinesin1.fasta
[KIFC1_Hs_Kinesin14]  length=673   saved → data/KIFC1_Hs_Kinesin14.fasta
```

### Step 3 — BLAST Hit Counts
```
Kinesin motor BLAST hits:  25–30
Kinesin-14 BLAST hits:     20–30
```
The exact number varies depending on PDB database updates at the time of the search.

### Step 6 — Merged Sequence Count
```
Fetched 15 sequences for KinesinMotor
Fetched 15 sequences for Kinesin14
Total unique sequences (both groups + queries): ~25–32
```

### Step 7 — Alignment Dimensions
```
Alignment: ~28 sequences × 800–1200 columns
```

### Step 8 — Trimming Result
```
Original length : ~1100 columns
Keeping columns : 12 → 980
Trimmed length  : ~968 columns
```

### Step 9 — Distance Matrix
- **Low distance (0.0–0.3)** between sequences within the same kinesin family
- **High distance (0.5–0.8)** between Kinesin-1 and Kinesin-14 sequences
- **Intermediate distances** between the two query sequences and their respective hits

### Step 10 — Deduplication
```
Input:  ~28 sequences
Output: ~24–28 sequences (removed 0–4 duplicates)
```

### Step 11 — ModelTest-NG Model Selection
```
Best model according to BIC: LG+I+G4

Model        LogL       AIC       BICc
LG+I+G4   -XXXXX    XXXXXX     XXXXXX   <- best
WAG+I+G4  -XXXXX    XXXXXX     XXXXXX
LG+G4     -XXXXX    XXXXXX     XXXXXX
```

### Step 12 — RAxML-NG Tree Inference
```
RAxML-NG v. 1.2.0
Analysis completed successfully
Final LogLikelihood: -XXXXX.XXXXXX
```

### Step 14 — Bootstrap Analysis
```
Bootstrap trees: 100
Bootstrap trees saved to: T2.raxml.bootstraps
```
Typically completes in 5–15 minutes on Colab CPU.

### Step 15 — Bootstrap Support Mapping
```
Bootstrap support mapped successfully → T3.raxml.support
```
- **>90%** at the node separating Kinesin-1 from Kinesin-14
- **>80%** within each family clade
- **50–80%** at shallower nodes within families

### Step 17 — Ancestral Reconstruction
```
Found: ASR.raxml.ancestralTree
Found: ASR.raxml.ancestralProbs
```
`ASR.raxml.ancestralProbs` is a tab-delimited matrix of marginal posterior probabilities across all 20 amino acids at each internal node.

---

## Key Biological Interpretation

The tree produced by this pipeline reproduces the core topology of the Duke Kinesin Resource tree. The two query proteins represent fundamentally different evolutionary solutions to microtubule-based motility:

| Feature | Kinesin Motor (Kinesin-1) | Kinesin-14 |
|---------|--------------------------|------------|
| Motor domain position | N-terminal | C-terminal |
| Directionality | Plus-end directed | Minus-end directed |
| Canonical members | KIF5A, KIF5B, KIF5C | Ncd, Kar3, KIFC1, Pkl1 |
| Primary function | Vesicle and organelle transport | Spindle assembly, chromosome segregation |
| Representative query | Human KIF5B (NP_004512) | Human KIFC1/HSET (NP_002254) |

The high bootstrap support at the node separating the two clades confirms they are monophyletic groups that diverged early in eukaryotic evolution.

---

## References

1. Kozlov, A. M., et al. (2019). RAxML-NG. *Bioinformatics*, 35(21), 4453–4455. https://doi.org/10.1093/bioinformatics/btz305
2. Katoh, K., and Standley, D. M. (2013). MAFFT Version 7. *Molecular Biology and Evolution*, 30(4), 772–780. https://doi.org/10.1093/molbev/mst010
3. Darriba, D., et al. (2020). ModelTest-NG. *Molecular Biology and Evolution*, 37(1), 291–294. https://doi.org/10.1093/molbev/msz189
4. Cock, P. J. A., et al. (2009). Biopython. *Bioinformatics*, 25(11), 1422–1423. https://doi.org/10.1093/bioinformatics/btp163
5. Duke Kinesin Resource. https://sites.duke.edu/kinesin/

---
