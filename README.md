# Richness
# Language: R
# Input: TXT
# Output: PREFIX
# Tested with: PluMA 1.1, R 4.0.0
# Dependencies: phyloseq 1.32.0, ape 5.4, psadd 0.1.2, ggplot2 3.3.1, microbiome 1.10.0

PluMA plugin that runs Principle Coordinate Analysis (Gower et al, 1966) to differentiate
sample sets.

The following are specified in the input TXT file, as tab-delimited keyword-value pairs.

otufile: OTU abundances (CSV)
mapping: Mapping table (CSV)
tree: Phylogenetic tree (CSV)
column: (STRING)
measure: (STRING) <- All will do all seven (Observed, Chao1, ACE, Shannon, Simpson, InvSimpson, Fisher)

The plot will then be generated using the provided prefix (prefix.pdf), with a supplemental prefix.csv file containing richness values.

Format of the otufile CSV:
- Rows are OTUs, columns are samples
- Entry (i,j) is the abundance of OTU in sample j.  Can be absolute or relative.

Format of the mapping CSV:
- Rows are samples, columns are metadata
- One single header row, i.e: sample,subject,visit

Format of the tree CSV;
- Rows are OTUs, columns are classifications
- One single header row, i.e. OTU,Kingdom,Phylum,Class,Order,Family,Genus,Species
