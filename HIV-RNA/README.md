Analysis Description
--------------------
Use a partition annotation file, and a time-scaled MCC tree, split the
alignment into partitions, fit partition-specific models with
substitution rates partitioned into terminal and interior nodes, report
the MLEs

- __Citation__: "_Combining chemical and evolutionary information to effectively evaluate HIV and SIV genetic heterogeneity and RNA structural conservation within and among infected hosts_", by Rife _et. al._

- __Written by__: Sergei L Kosakovsky Pond

- __Contact Information__: spond@temple.edu

- __Analysis Version__: 1.00

## How to run

1. Download and install HyPhy (www.hyphy.org). This script is confirmed to work with v2.3.3 (https://github.com/veg/hyphy/releases/tag/v2.3.3), but may also work with newer versions
2. Run `fit-nexus-partitions.bf` in HyPhy, e.g. <code>$HYPHYMP /path/to/fit-nexus-partitions.bf</code>
3. When prompted, 
   - supply the name of the NEXUS file with partitions specified as CHARSET blocks (see `data/N01.nex` for an example)
   - Agree to use the tree in the file or specify a separate tree (in Newick format with branch lengths in units of time)
4. Output is written to standard out as a Markdown text

## Example output 


### Loaded an MSA with **443** sequences and **1308** sites from `/Volumes/sergei-raid/Coding/pubs/HIV-RNA/hyphy-2.3.3/../data/N01.nex`
Read **10** partitions on **1308** sites.

### Fitting partition CODON1+2 with **358** constant sites, **94** variable sites, and **47** phylogenetically informative sites
* log (L)  **-3476.108**
* &pi; (A)  **0.232**
* &pi; (C)  **0.228**
* &pi; (G)  **0.292**
* &pi; (T)  **0.248**
* Internal node rate  **0.0000176 [0.0000150 - 0.0000208]**
* Internal node kappa **5.863 [4.710 - 7.210]**

### Fitting partition CODON3 with **126** constant sites, **111** variable sites, and **58** phylogenetically informative sites
* log (L)  **-3078.507**
* &pi; (A)  **0.212**
* &pi; (C)  **0.195**
* &pi; (G)  **0.266**
* &pi; (T)  **0.326**
* Internal node rate  **0.0000330 [0.0000272 - 0.0000390]**
* Internal node kappa **7.140 [5.746 - 8.747]**

### Fitting partition CODON4+5 with **310** constant sites, **108** variable sites, and **48** phylogenetically informative sites
* log (L)  **-2952.911**
* &pi; (A)  **0.538**
* &pi; (C)  **0.107**
* &pi; (G)  **0.161**
* &pi; (T)  **0.194**
* Internal node rate  **0.0000138 [0.0000114 - 0.0000175]**
* Internal node kappa **4.997 [3.846 - 6.381]**

### Fitting partition CODON6 with **116** constant sites, **79** variable sites, and **47** phylogenetically informative sites
* log (L)  **-2602.324**
* &pi; (A)  **0.384**
* &pi; (C)  **0.111**
* &pi; (G)  **0.130**
* &pi; (T)  **0.375**
* Internal node rate  **0.0000335 [0.0000289 - 0.0000383]**
* Internal node kappa **32.190 [26.282 - 40.115]**