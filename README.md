# Comparative gene regulatory network mapping of Brassicaceae members with differential drought tolerance

This repository contains the analytical codebase, visualization pipelines and network modeling scripts developed to analyze the cross-species drought-responsive gene regulatory networks (GRNs) within the Brassicaceae family.

## File registry & Pipeline mapping
1. Network topology

`Bipartite Package Nestedness and Connectedness/bipartite_indices.r`
- Translates raw cross-species Y1H interactions into directed bipartite network graphs.
- Network connectance, niche overlap and modularity indices were calculated to evaluate the shared and unique regulatory spaces.

`Bipartite Package Nestedness and Connectedness/degree_distribution.r`
- Establishes un-directed bipartite network models to compute node degrees for both individual TFs and downstream target promoters across all three species.

2. Coregulation & Expression evolutionary analysis

`Coexpression_Analysis/linear_correlation.py`
- Integrates time-course expression profiles with physical protein-DNA interactions to determine if interacting TF-promoter pairs exhibit tighter coregulatory linear relationships than non-interacting background pairs.

3. Genome-Wide promoter parsing & Motif architecture

`Genome Motif analysis/Genome-wide-promoter-parsing.ipynb`
- A notebook that retrieves, processes and normalizes whole-genome sequence assemblies and structural GFF3 annotations directly from Ensembl Plants and NCBI repositories.
