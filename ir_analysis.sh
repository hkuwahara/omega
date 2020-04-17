#! /bin/bash


sample_dir=$1
chr=$2
strand=$3
read_len=$4
annotation_dir=$5

id="${chr}${strand}"

output_dir=${sample_dir}/sj_exon_map
mkdir -p ${output_dir}

BAMFILE=${sample_dir}/bam/${id}.bam
NONSPLIT_BAMFILE=${sample_dir}/bam/nonsplit-${id}.bam

samtools view -h ${BAMFILE} | awk 'BEGIN{FS="\t"; OFS="\t"} $1 ~/@/{print $0; next;}  !($6 ~/N/){print $0; next;}' | samtools sort | samtools view -Sb > ${NONSPLIT_BAMFILE}
bedtools coverage -a ${annotation_dir}/nonexonic-${id}.bed -b ${NONSPLIT_BAMFILE} | awk -v readLen="${read_len}" 'BEGIN{FS="\t"; OFS="\t";} {print $1, $2, ($3+1), $4, $5, $6, (readLen * $4 / $6); next;}' > ${output_dir}/nonexonic-${id}.cov


exon_exon_pair=${annotation_dir}/coding_exon_exon_juncs_${id}.tsv
ir_count=${output_dir}/ir-count-${id}.tsv

if [ "${strand}" = "+" ]; then	
	nonexonic_cov=${SJ_EXON_MAP_DIR}/nonexonic-${id}.cov
	awk 'BEGIN{FS="\t"; OFS="\t"; print "donor-site", "acceptor-site", "ir-count";} NR==FNR{ a[$2 OFS $3] = $7; next;} {if(!(($5 OFS $6) in b)){x=a[$3 OFS $4] + 0; print $3, $4, $5, $6, x; b[$5 OFS $6] = 1;}}' ${output_dir}/nonexonic-${id}.cov ${exon_exon_pair} > ${ir_count} 
else	
	awk 'BEGIN{FS="\t"; OFS="\t"; print "donor-site", "acceptor-site", "ir-count";} NR==FNR{ a[$3 OFS $2] = $7; next;} {if(!(($5 OFS $6) in b)){x=a[$3 OFS $4] + 0; print $3, $4, $5, $6, x; b[$5 OFS $6] = 1;}}' ${output_dir}/nonexonic-${id}.cov ${exon_exon_pair} > ${ir_count} 
fi

rm -f ${NONSPLIT_BAMFILE}
rm -f ${output_dir}/nonexonic-${id}.*

