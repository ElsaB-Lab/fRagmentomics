# Rule for adding the mutated or wild-type status of fragments in different samples.
rule add_mutations_status:
    input:
        tsv='%s/fragments/fragmentomics/{sample}_{chr}_{pos}.tsv' % R_FOLDER,
    output:
        tsv='%s/fragments/fragmentomics_mut_status/{sample}_{chr}_{pos}_mut_status.tsv' % R_FOLDER,
    benchmark:
        "%s/fragments/fragmentomics_mut_status/{sample}_{chr}_{pos}_mut_status.tsv" % B_FOLDER
    log:
        "%s/fragments/fragmentomics_mut_status/{sample}_{chr}_{pos}_mut_status.log" % L_FOLDER
    conda:
        "../envs/r.yaml"
    threads: 1
    resources:
        queue="shortq",
        mem_mb=500,
        time_min=20
    shell:
        """
        Rscript workflow/scripts/04.1_mutations_status.R {input.tsv} {wildcards.pos} {output.tsv}
        """
