#! /bin/bash


sample_dir=$1
chr=$2
strand=$3
annotation_dir=$4

id="${chr}${strand}"

output_dir=${sample_dir}/sj

mkdir -p ${output_dir}

BAMFILE="${sample_dir}/bam/${id}.bam"
SJ_FILE="${output_dir}/${id}"

awk 'BEGIN{FS="\t"; OFS="\t";} NR==1{ch=$1;} !($4 in lexon){lexon[$4] = $2; rexon[$4] = $3; next;} {ljunc = rexon[$4] + 1; rjunc = $2; juncs[ljunc OFS rjunc] += 1; lexon[$4] = $2; rexon[$4] = $3; 		next;} END{for( j in juncs ){ count = juncs[j]; if( count >= 5 ) {print ch, j, ".", ".", ".", count, ".", ".";} } }' <(samtools view -h ${BAMFILE} | awk 'BEGIN{FS="\t"; OFS="\t";} ($1 ~/@/){print $0; next;} ($6 ~ /N/) && ($5 == 255) {print $0}' | samtools view -bS - | bamToBed -bed12 | bed12ToBed6 -i stdin | sort -k2,3n) | sort -k2,3n > ${SJ_FILE}


Rscript --vanilla ./map_junctions_to_exon_pairs.R ${chr} ${strand} ${sample_dir} ${annotation_dir}
