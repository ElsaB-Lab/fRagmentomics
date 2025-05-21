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

##### Helper functions #####

def get_data_irods(wildcards):
    """Get irods paths for data of given sample."""
    table_sam = sam_table.loc[wildcards.sample, ["BAM_Path", "BAI_Path", "PDF_Path", "XML_Path"]].dropna()
    return {"BAM": table_sam.BAM_Path,
            "BAI": table_sam.BAI_Path,
            "PDF": table_sam.PDF_Path,
            "XML": table_sam.XML_Path}
