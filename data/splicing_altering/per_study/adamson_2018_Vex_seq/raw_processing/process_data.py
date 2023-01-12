import pandas as pd

target_variants=pd.read_csv('all_variants.txt', sep="\t")
dpsi_df = pd.read_csv('dPSI_values.txt', sep="\t")

# Select variants with dPSI less than 1%
dpsi_no_change = dpsi_df[dpsi_df.delta_PSI.abs() < 1]
dpsi_no_change = pd.merge(dpsi_no_change, target_variants, on='variant', how='left')

# Select variants with dPSI higher than 10%
#dpsi_change =  dpsi_df[dpsi_df.delta_PSI.abs() > 20]
#dpsi_change = pd.merge(dpsi_change, target_variants, on='variant', how='left')

# Remove variants in UTRs and non-coding transcripts
dpsi_no_change = dpsi_no_change[~dpsi_no_change.annotation.isin(['3\' UTR', '5\' UTR', 'non coding transcript exon'])]
dpsi_no_change = dpsi_no_change[dpsi_no_change.variant != "chr16_67690790_G_A"]

# Split exonic and intronic
exonic_no_change = dpsi_no_change[dpsi_no_change.annotation_no_SR != "intron"]
intronic_no_change = dpsi_no_change[dpsi_no_change.annotation_no_SR == "intron"]
#intronic_change = dpsi_change[dpsi_change.annotation_no_SR == "intron"]

# Create VEP input
for i, df in enumerate([exonic_no_change, intronic_no_change]): #, intronic_change]):
    if i == 0:
        outfile="1_exonic_to_vep.tsv"
    elif i == 1:
        outfile="1_intronic_no_change_to_vep.tsv"
#    elif i == 2:
#        outfile="1_intronic_with_change_to_vep.tsv"
    out = open(outfile, 'w')
    for _, row in df.iterrows():
        v_id = row.variant.split("_")
        chrom=v_id[0].replace("chr","")
        pos=v_id[1]
        ref=v_id[2]
        alt=v_id[3]
        if len(ref) == len(alt):
            out.write('{}\t{}\t{}\t{}/{}\t+\t{}\n'.format(chrom, pos, pos, ref, alt, row.variant))
        elif len(ref) > len(alt):
            diff = len(ref) - len(alt) 
            start_pos = int(pos) + 1
            end_pos = int(pos) + diff
            out_str = '{}\t{}\t{}\t{}/{}\t+\t{}\n'.format(chrom, start_pos, end_pos, ref[1:], "-", row.variant)
            out.write(out_str)
        else:
            print(v_id)

