#! /bin/bash


sample_dir=$1
chr=$2
strand=$3
annotation_dir=$4

id="${chr}${strand}"

exon_tran_map=${annotation_dir}/coding_exon_trans_map_${id}.tsv
sj_exon_pair_tran_map=${annotation_dir}/sj_exon_pair_trans_${id}.tsv

sj_file=${sample_dir}/sj_exon_map/${id}.tsv
ir_file=${sample_dir}/sj_exon_map/ir-count-${id}.tsv

normal_junc_fraction=${sample_dir}/sj_exon_map/normal_junc_fraction_${id}.tsv
fraction_normal_tran_map=${sample_dir}/sj_exon_map/frac_normal_trans_${id}.tsv

quant_data=${sample_dir}/quant/transcript_abundance_${id}.tsv
adjusted_quant_data=${sample_dir}/quant/adjusted_quant_${id}.tsv


if [ "${strand}" = "+" ]; then	
	awk 'BEGIN{FS="\t"; OFS="\t";} FNR==1{next;} NR==FNR{abnormal[$3 OFS $4] = $5; next;} {a[$4 OFS $6] = 1; if($5 != 0 || $7 != 0){abnormal[$4 OFS $6] += $3;} else {normal[$4 OFS $6] += $3;} next;} END{for(key in a){print normal[key]/(abnormal[key] + normal[key]), key}}' ${ir_file} ${sj_file} > ${normal_junc_fraction}  
else
	awk 'BEGIN{FS="\t"; OFS="\t";} FNR==1{next;} NR==FNR{abnormal[$3 OFS $4] = $5; next;} {a[$6 OFS $4] = 1; if($5 != 0 || $7 != 0){abnormal[$6 OFS $4] += $3;} else {normal[$6 OFS $4] += $3;} next;} END{for(key in a){print normal[key]/(abnormal[key] + normal[key]), key}}' ${ir_file} ${sj_file} > ${normal_junc_fraction}  
fi

awk 'BEGIN{FS="\t"; OFS="\t"; frac = 1;} NR==FNR{a[$2 OFS $3] = $1; next;} (($1 OFS $2 OFS $3) in b){w = (($6 OFS $7) in a) ? a[$6 OFS $7] : 1; b[$1 OFS $2 OFS $3] = b[$1 OFS $2 OFS $3] * w; next;} {b[$1 OFS $2 OFS $3] = (($6 OFS $7) in a) ? a[$6 OFS $7] : 1;} END{for( key in b ){ print key, b[key];}}'  ${normal_junc_fraction} ${sj_exon_pair_tran_map} > ${fraction_normal_tran_map}

awk 'BEGIN{FS="\t"; OFS="\t"; print "start", "end", "ensembl_gene_id", "ensembl_transcript_id", "transcript_symbol", "gene_symbol", "length", "estimated_length", "count", "XXX", "normal_frac", "adjusted_count";} NR==FNR{ a[$2] = $4; next;} FNR>1{ adjusted = ( $4 in a ) ? a[$4] OFS (a[$4] * $9) : "1" OFS $9; print $0, adjusted;}'  ${fraction_normal_tran_map} ${quant_data} | cut --complement -f10 > ${adjusted_quant_data}


