# CNomplexity
Package for investigating complex copy number changes and structural variants

Dependencies: parallel, GenomicRanges

## For identifying putative chromothriptic regions: chromothripsis()
Implements three of the tests for chromothripsis from [Korbel & Campbell (2013)](http://www.sciencedirect.com/science/article/pii/S0092867413002122). 
Tests for departure from exponentially distributed distances between breakpoints, equal proportion of HH, HT, TH and TT fusions partners, and random choice of fusion partners.

## For identifying genome doubled samples: genomeDoubling()
Implements tests for genome doubling from [Dewhurst et al. (2013)](http://cancerdiscovery.aacrjournals.org/content/4/2/175.short) and [Carter et al. (2012)](https://www.nature.com/nbt/journal/v30/n5/full/nbt.2203.html?foxtrotcallback=true).

## For scoring chromosome level copy number changes: CNscore()
Implements copy number scoring from [Davoli et al. (2017)](http://science.sciencemag.org/content/355/6322/eaaf8399).
