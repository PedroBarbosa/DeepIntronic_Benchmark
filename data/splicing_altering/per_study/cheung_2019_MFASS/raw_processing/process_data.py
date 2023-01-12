import pandas as pd

variants=pd.read_csv('variants.txt', sep="\t")

# Select variants with dPSI less than 2% in both replicates
dpsi_df = variants[(variants.v1_dpsi.abs() < 0.02) & (variants.v2_dpsi.abs() < 0.02)]

# Select variants with dPSI larger than 20% in both replicates
#dpsi_with_change = variants[(variants.v1_dpsi.abs() >= 0.1) & (variants.v2_dpsi.abs() >= 0.1)]

# Downstream intronic 
downstream = dpsi_df[dpsi_df.label == "downstr_intron"]

# Upstream PY tract 
upstream_py = dpsi_df[(dpsi_df.label == "upstr_intron") & ((dpsi_df.rel_position_feature < -44) | (dpsi_df.rel_position_feature > -18) & (dpsi_df.rel_position_feature < -2))] 
upstream_branchpoint = dpsi_df[(dpsi_df.label == "upstr_intron") & (dpsi_df.rel_position_feature < -18) & (dpsi_df.rel_position_feature >= -44)]

# Create VEP input
all_exon_ids = []
for i, df in enumerate([downstream, upstream_py, upstream_branchpoint]):
    if i == 0:
        outfile="1_downstream_to_vep.tsv"
    elif i == 1:
        outfile="1_upstream_py_tract_to_vep.tsv"
    elif i == 2:
        outfile="1_upstream_branchpoint_to_vep.tsv"

    out = open(outfile, 'w')
    for _, row in df.iterrows():
        chrom=row.chr.replace("chr","")
        pos=str(row.snp_position)
        ref=row.ref_allele
        alt=row.alt_allele
        strand=row.strand
        label=chrom + "_" + pos + "_" + ref + "_" + alt + "_" + row.ensembl_id
        all_exon_ids.append(row.ensembl_id)
        out.write('{}\t{}\t{}\t{}/{}\t{}\t{}\n'.format(chrom, pos, pos, ref, alt, "+", label))

out_exon_ids=open('1_exon_ids.txt', 'w')
[out_exon_ids.write(x + "\n") for x in set(all_exon_ids)]
