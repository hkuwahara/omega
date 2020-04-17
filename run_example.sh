#! /bin/bash

toplevel=./example
chr="chr21"
strand="-"

##
## extract sj list from bam file and map sjs to closest annotated junctions.  
## the splicing junction file is located at subdirectory sj, while the mapping analysis file is located at sj_exon_map
##
./sj_analysis.sh ${toplevel} ${chr} ${strand} ${toplevel}/annotations 

##
## from non-split reads, compute the weight for intron retention events and generate a list under sj_exon_map
##
./ir_analysis.sh ${toplevel} ${chr} ${strand} 150 ${toplevel}/annotations 



#
# adjust the transcript-level quantification using the weight computed via sj and ir analysis.  
#
./adjust_quant.sh ${toplevel} ${chr} ${strand} ${toplevel}/annotations 

#
# to get normalized quant in CPM, first compute the scaling factor.  this should have been computed before running the main script via command like below:
#
#scaling=$(awk 'BEGIN{FS="\t";} FNR>1{sum += $9; next;} END{print 1000000/sum;}' ${toplevel}/quant/transcript_abundance_chr*.tsv)
scaling="0.0134785"
 
#
# finally, compute cpm, omega, and rate.
#
./omega_gen.sh ${toplevel} ${chr} ${strand} ${scaling} 


