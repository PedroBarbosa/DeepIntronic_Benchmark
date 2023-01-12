from numpy import var
import pandas as pd


def assign_mut_class(row: pd.Series):

    if row['cryptic SSs'] == 'donor':

        if row['location relative to cryptic SSs'] == 'downstream':

            if row['distance to cryptic SSs'] > 3:
                row['mutation_class']  = 'strenghtening_donor'

            else:
                row['mutation_class']  = 'new_splice_donor'
        
        elif row['location relative to cryptic SSs'] == 'upstream':
            if row['distance to cryptic SSs'] > 2:
                row['mutation_class']  = 'change_sre'
            else:
                row['mutation_class'] = 'new_splice_donor'
        else:
            print ('Error: {}'.format(row['location relative to cryptic SSs'] ))
            exit(1)
    
    elif row['cryptic SSs'] == 'acceptor':
        if row['location relative to cryptic SSs'] == 'upstream':

            if row['distance to cryptic SSs'] < 3:
                row['mutation_class']  = 'new_splice_acceptor'             

            elif 3 < row['distance to cryptic SSs'] < 18:
                row['mutation_class']  = 'strenghtening_acceptor'

            elif 18 <= row['distance to cryptic SSs'] <= 44:
                row['mutation_class']  = 'branchpoint_associated'
            else:
                print('Error. Distant mutation: {}'.format(row['distance to cryptic SSs']))
                exit(1)

        elif row['location relative to cryptic SSs'] == 'downstream':
            if row['distance to cryptic SSs'] > 2:
                row['mutation_class'] = 'change_sre'
            else:
                row['mutation_class']  = 'new_splice_acceptor'
        
        else:
            print('Error: {}'.format(row['location relative to cryptic SSs'] ))
            exit(1)

    return row

variants_to_eval = pd.read_csv('2_deep_tabular_after_VEP.tsv', sep="\t")
original_to_extract_class = pd.read_csv('0_deep_intronic_original.tsv', sep='\t')[['chr', 'start', 'strand', 'cryptic SSs', 'location relative to cryptic SSs', 'distance to cryptic SSs', 'functional_consequence' ]]

original_to_extract_class = original_to_extract_class.rename(columns={'chr': '#CHROM', 'start': 'POS'})

joined = pd.merge(variants_to_eval, original_to_extract_class, on=['#CHROM', 'POS'], how='left')
joined = joined.apply(assign_mut_class, axis=1, result_type='expand')
joined['Source'] = 'jung_2021'

cols = ['#CHROM', 'POS', 'REF', 'ALT', 'ID', 'Gene', 'SYMBOL', 'HGVSc', 'HGVSg', 'INTRON', 'gnomADg_AF', 'mutation_class', 'Source']
joined[cols].to_csv('deep_tabular_final.tsv', sep="\t", index=False)
