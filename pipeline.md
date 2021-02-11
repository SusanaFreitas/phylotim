This rep will keep all info about the phylogenetic analyses I will do for Timema.

##########################################################################################################################3

Pipeline
0 - get reads (OMG, this was a bit messy)
1 - filter reads
2 - map reads against Tce assembly


###################
## 0 - get reads ##

Reads were from mixed sources. Some reads were RADseq done on Chloe samples, others were Nosil GBS or whole genome sequencing.
GBS samples from the Nosil group were more or less 3x less (in terms total number of basis) than the Schwander group. 
(from a quick and dirty estimate)
So I will pool the Nosil GBS samples in groups of 3.

Nosil reads were trimmed to 60 bp
Schwander reads trimmed to 80 bp

"""bash
cat inds_trim80 | while read line; do co reads/"$line"* . ; done
"""
