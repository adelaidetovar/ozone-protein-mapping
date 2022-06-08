## A locus on chromosome 15 contributes to acute ozone-induced lung injury in Collaborative Cross mice

This repository contains analyses performed for the accompanying article, and were originally written and performed using R 4.1.1. The analysis is separated into several files, all sourcing a file 
with the 
environment and accessory functions. Each should be run sequentially. A separate folder contains necessary supplemental data files.

To recreate analyses, clone this repo and work from the analysis directory. Before performing QTL mapping, you will need to create a genome cache using the code within `analysis/3a_genome_cache.R` -- this is sourced when running the mapping code. Retrieve [Collaborative Cross MRCAs](http://csbio.unc.edu/CCstatus/index.py?run=FounderProbs) for all strains included in the study and [MegaMUGA reference clusters](http://csbio.unc.edu/CCstatus/Media/snps.megamuga.Rdata) and adjust paths in code accordingly. 

This analysis is open source under an MIT license.
