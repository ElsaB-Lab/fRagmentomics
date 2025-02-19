import argparse
import numpy as np
import pandas as pd
import re
import functools
print = functools.partial(print, flush=True)

def print_status_changes(df_table):
    old_ok = ~df_table["IRODS_Status_Old"].isnull()
    new_ok = ~df_table["IRODS_Status"].isnull()
    total = df_table.shape[0]
    changed_to_ok = ~old_ok & new_ok
    changed_to_nok = old_ok & ~new_ok
    print("%d/%d samples changed from Not OK to OK:" % (sum(changed_to_ok), total))
    if sum(changed_to_ok)>0:
        sample_ids_changed_to_ok = df_table.loc[changed_to_ok, "Sample_Id"].tolist()
        print("\t" + "\n\t".join(sample_ids_changed_to_ok))
    print("%d/%d samples changed from OK to Not OK" % (sum(changed_to_nok), total))
    if sum(changed_to_nok)>0:
        sample_ids_changed_to_nok = df_table.loc[changed_to_nok, "Sample_Id"].tolist()
        print("\t" + "\n\t".join(sample_ids_changed_to_nok))
    return sum(changed_to_ok), sum(changed_to_nok)


def main(args):
    df_table = pd.read_table(args.table)
    if "IRODS_Status" not in df_table:
        raise ValueError("cannot check IRODS status as 'IRODS_Status' column is absent from %s" % args.table)
    else:
        df_table = df_table.rename(columns={"IRODS_Status": "IRODS_Status_Old"})

    df_irods = pd.read_table(args.irods, header=None)
    df_irods.columns = ["Dataset"]
    mask_keep = df_irods["Dataset"].apply(lambda x: x.startswith("collection:"))
    df_irods = df_irods.loc[mask_keep].copy()
    df_irods["Dataset"] = df_irods["Dataset"].apply(lambda x: re.sub("collection: ", "", x))
    df_irods["IRODS_Status"] = "OK"

    col_x = "BAM_Path"
    col_y = "Dataset"
    df_table[col_y] = df_table[col_x].str.extract(r'^(.+)(?=\/archive)')
    df_table = df_table.merge(df_irods, how="left", on="Dataset")

    n_changed_to_ok, n_changed_to_nok = print_status_changes(df_table)

    if n_changed_to_ok + n_changed_to_nok > 0:
        raise ValueError("some files changed status after automatic IRODS check. Please review your " + \
                         "the values of 'IRODS_Status' in your config file 'samples.tsv' or disable automatic " + \
                         "IRODS check in the config.yaml file.")
    else:
        open(args.output, 'w').close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Checks whether all files dataset are still
                                     available in IRODS or not""")
    parser.add_argument("--table", type=str, help="Path to table of samples")
    parser.add_argument("--irods", type=str, help="Path to output of imeta qu -C command.")
    parser.add_argument("--output", type=str, help="Path to output table with status.")

    args = parser.parse_args()
    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n", end="")

    main(args)
