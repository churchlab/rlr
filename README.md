# rlr
Retron Library Recombineering
Maintained by Max Schubert mgschubert@gmail.com and Daniel Goodman dbgoodman@gmail.com .

This repository contains all scripts used to generate results for this publication, and the processed data required to run these scripts.

For raw sequencing Data, see the submission to the NCBI GEO database linked to this publication.

The Majority of figures are generated in R, and rMarkdown(.rmd) files are given for each figure or figure panel. In some cases one rMarkdown file is used to generate several figure panels. Scripts within these files retreive data from folders within this repository, generally organized by the figure to which they pertain.

Population simulations and simulations of editing are performed in Matlab, giving rise to images in supplemental figure 2 and figure 2E. These scripts can be found in supplemental_figures/fig_s2

Alignment and calculation of coverage depth in figures 5B, 5C, and S6B, as well as deriving the "template lengths" depicted in S6A, are performed as described in the Online Methods section. Clarification and example scripting commands are available upon request.

