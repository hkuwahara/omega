#! /bin/bash


sample_dir=$1
chr=$2
strand=$3
scaling=$4

id="${chr}${strand}"
cpm_dir=${sample_dir}/cpm
omega_dir=${sample_dir}/omega

mkdir -p ${cpm_dir}
mkdir -p ${omega_dir}

quant_data=${sample_dir}/quant/transcript_abundance_${id}.tsv
adjusted_quant_data=${sample_dir}/quant/adjusted_quant_${id}.tsv

cpm_data=${cpm_dir}/transcript_abundance_${id}.tmp
adjusted_cpm_data=${omega_dir}/adjusted_cpm_${id}.tmp

gene_cpm_data=${cpm_dir}/cpm_${id}.tsv
omega_data=${omega_dir}/omega_${id}.tsv
rate_data=${omega_dir}/rate_${id}.tsv


printf "chr\tstart\tend\tstrand\tensembl_gene_id\tensembl_transcript_id\ttranscript_symbol\tgene_symbol\tlength\testimated_length\tcount\tCPM\n" > ${cpm_data}
awk -v ch="${chr#chr}" -v st="${strand}" -v scale_factor="${scaling}" 'BEGIN{FS="\t"; OFS="\t";} NR==1{next;} {val=$9 * scale_factor; print ch, $1, $2, st, $3, $4, $5, $6, $7, $8, $9, val;}' ${quant_data} >>  ${cpm_data}  

printf "chr\tstart\tend\tstrand\tensembl_gene_id\tensembl_transcript_id\ttranscript_symbol\tgene_symbol\tlength\testimated_length\tadjusted_count\tadjusted_CPM\n" > ${adjusted_cpm_data}
awk -v ch="${chr#chr}" -v st="${strand}" -v scale_factor="${scaling}" 'BEGIN{FS="\t"; OFS="\t";} NR==1{next;} {val=$11 * scale_factor; print ch, $1, $2, st, $3, $4, $5, $6, $7, $8, $11, val;}' ${adjusted_quant_data} >>  ${adjusted_cpm_data}  


awk 'BEGIN{FS="\t"; OFS="\t"; print "chr", "strand", "ensembl_gene_id", "gene_symbol", "CPM";} NR==1{next;} {a[$5] += $12; if(!($5 in b)){b[$5] = $1 OFS $4 OFS $5 OFS $8;}} END{for(g in a){print b[g], a[g];}}' ${cpm_data} | sort -k3,3 > ${gene_cpm_data} 


awk 'BEGIN{FS="\t"; OFS="\t"; print "chr", "strand", "ensembl_gene_id", "gene_symbol", "omega";} NR==1{next;} {a[$5] += $12; if(!($5 in b)){b[$5] = $1 OFS $4 OFS $5 OFS $8;}} END{for(g in a){print b[g], a[g];}}' ${adjusted_cpm_data} | sort -k3,3 > ${omega_data} 

awk 'BEGIN{FS="\t"; OFS="\t"; print "chr", "strand", "ensembl_gene_id", "gene_symbol", "rate";} FNR==1{next;} NR==FNR{a[$3] = $5; next;} {x = a[$3]; print $1, $2, $3, $4, (x > 1e-6) ? $5/x : 1;}' ${gene_cpm_data} ${omega_data} > ${rate_data}  

rm -f ${cpm_dir}/*.tmp
rm -f ${omega_dir}/*.tmp ${adjusted_quant_data}

