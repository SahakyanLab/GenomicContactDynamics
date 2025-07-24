## Analysis of long-range contacts across cell types outlines a core sequence determinant of 3D genome organisation
Liezel Tamon, Zahra Fahmi, James Ashford, Rosana Collepardo-Guevara, and Aleksandr B. Sahakyan*

Each directory contains a README explaining each script. Directories relevant for each main figure:

#### Figure 1. Contact persistence to isolate and investigate the core genomic contacts, their determinants and implications.
- C. `11_Complementarity` to generate contact maps.

#### Figure 2. Visualisation of persistent contacts and the chromosome organisation they mediate.
- A,C. `4_ArcPlot` to generate arc plots representing contacts per chromosome.
- B,D. [CoreGenomeExplorer Shiny app](https://github.com/liezeltamon/CoreGenomeExplorerLite) to generate network plots representing contacts per chromosome.

#### Figure 3. Persistent contacts enriched for contacts with features associated with heterochromatin and preferential AT sequence composition.
- A. `8_FeatureVsPersist/A5_metaplot.R`
- B. `7_FeaturePermutation/D1_foivsij.R`
- C. `23_Hub_Gene_Expression/E1_geneExprVsCp.R`
- D. `16_GeneVsPersist`
- E. `16_GeneVsPersist`
- F. `20_ReplicationTiming`
- G. `19_MutationRatesVsPersist`
  
#### Figure 4. Higher sequence complementarity between persistent contacts.
- `11_Complementarity` to calculate sequence complementarity and generate contact maps.
- `12_Shuffling` to generate shuffled contacts.
  
#### Figure 5. Genome-wide cII values recapitulate some Hi-C features in additional human datasets of different resolutions and in another species, Drosophila.
- `11_Complementarity` to calculate sequence complementarity and generate contact maps.

#### Figure 6. Repeat contribution to the observed sequence complementarity of contacts.
- `18_RepeatVsPersist`

#### Others
- `lib` contains utility functions used in scripts across directories
