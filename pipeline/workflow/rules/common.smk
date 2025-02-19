from itertools import product
import pandas as pd
import re
from snakemake.utils import validate
from snakemake.utils import min_version

min_version("5.4.0")

B_FOLDER = "workflow/benchmarks"
L_FOLDER = "workflow/logs"
R_FOLDER = "results"

###### Config file and sample sheets #####
configfile: "config/config.yaml"

# samples table
sam_table = pd.read_table(config["general"]["samples"], dtype=str).set_index(["Sample_Id"], drop=False)
if "IRODS_Status" in sam_table:
    mask_ok = sam_table["IRODS_Status"]=="OK"
else:
    mask_ok = pd.Series(True, index=sam_table.index)
samples = sam_table.loc[mask_ok]["Sample_Id"].tolist()


# mutations table
mut_table = pd.read_table(config["general"]["mutations"])
chr_all = mut_table["Chromosome"].unique().tolist()
pos_all = mut_table["Position"].unique().tolist()
#alt_all = mut_table["Alt"].unique().tolist()

##### Helper functions #####

def filter_combinator(combinator, comblist, white_list=True):
    def filtered_combinator(*args, **kwargs):
        for wc_comb in combinator(*args, **kwargs):
        # Use frozenset instead of tuple
        # in order to accomodate
        # unpredictable wildcard order
            if white_list:
                if frozenset(wc_comb) in comblist:
                    yield wc_comb
            else:
                if frozenset(wc_comb) not in comblist:
                    yield wc_comb
    return filtered_combinator

def get_allowed_sample_chr_pos():
    allowed = []
    df_mut = pd.read_table(config["general"]["mutations"])
    for (sample, chr, pos) in zip(df_mut["Sample_Id"], df_mut["Chromosome"], df_mut["Position"]):
        allowed.append(frozenset({("sample", sample), ("chr", chr), ("pos", pos)}))
    return filter_combinator(product, allowed, white_list=True)

#def get_allowed_sample_chr_pos():
#    allowed = []
#    df_mut = pd.read_table(config["general"]["mutations"])
#    for (sample, chr, pos, alt) in zip(df_mut["Sample_Id"], df_mut["Chromosome"], df_mut["Position"], df_mut["Alt"]):
#        allowed.append(frozenset({("sample", sample), ("chr", chr), ("pos", pos), ("alt", alt)}))
#    return filter_combinator(product, allowed, white_list=True)


def get_data_irods(wildcards):
    """Get irods paths for data of given sample."""
    table_sam = sam_table.loc[wildcards.sample, ["BAM_Path", "BAI_Path", "PDF_Path", "XML_Path"]].dropna()
    return {"BAM": table_sam.BAM_Path,
            "BAI": table_sam.BAI_Path,
            "PDF": table_sam.PDF_Path,
            "XML": table_sam.XML_Path}

# Functions to add wildcard to call fragmentomics
def get_ref(sample, chr, pos):
    pos = int(pos)
    row = mut_table.query("Sample_Id == @sample and Chromosome == @chr and Position == @pos")
    if len(row) == 1:
        return row["Ref"].values[0]
    else:
        return "NA"

def get_alt(sample, chr, pos):
    pos = int(pos)
    row = mut_table.query("Sample_Id == @sample and Chromosome == @chr and Position == @pos")
    if len(row) == 1:
        return row["Alt"].values[0]
    else:
        return "NA"

def get_mutation_type(sample, chr, pos):
    pos = int(pos)
    row = mut_table.query("Sample_Id == @sample and Chromosome == @chr and Position == @pos")
    if len(row) == 1:
        return row["Mutation_type"].values[0]
    else:
        return "NA"

def get_del_info(sample, chr, pos):
    pos = int(pos)
    row = mut_table.query("Sample_Id == @sample and Chromosome == @chr and Position == @pos")
    if len(row) == 1:
        return row["Del_info"].values[0]
    else:
        return "NA"
