#!/usr/bin/env python

import pandas as pd
from collections import defaultdict
import os

sample = "$sample"
hits_file = "$hits"
results_file = sample + "tax_counts.txt"

lineage_levels = ['Superkingdom', 'Phylum', 'Class', 'Order',
                  'Family', 'Genus', 'Species']

level_index_map = {lvl: i for i, lvl in enumerate(lineage_levels)}

def extract_lineage_level(lineage_str, level_idx):
    """Efficiently extract lineage level by index."""
    if level_idx == 0:
        return lineage_str.replace(" ","").split(";")[level_idx]
    
    else:
        return ";".join(lineage_str.replace(" ","").split(";")[0:level_idx + 1])

def process_file(hits_file, chunksize = 1000000):
    lldict = {lvl: defaultdict(int) for lvl in lineage_levels}

    for chunk in pd.read_csv(hits_file, compression='gzip', sep="\t", chunksize=chunksize, 
                             skiprows=0, names = header):
        
        mappedchunk = chunk[chunk['status'] == 'C']

        for level in lineage_levels:
            idx = level_index_map[level]
            values = mappedchunk['taxonomic_lineage'].apply(lambda x: extract_lineage_level(x, idx))
            counts = values.value_counts()

            ldict = lldict[level]

            for taxon, count in counts.items():
                ldict[taxon] += count
    return lldict

def lldict_to_dataframe(lldict, sample_id):
    records = []
    for level, taxon_dict in lldict.items():
        for taxon, count in taxon_dict.items():
            records.append((sample_id, level, taxon, count))
    return pd.DataFrame(records, columns=['sample_id', 'lineage_level', 'taxon', 'count'])


lldict = process_file(hits_file)
df = lldict_to_dataframe(lldict, sample)
df.to_csv(file = results_file, sep = "\t")
