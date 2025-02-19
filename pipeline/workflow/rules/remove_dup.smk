# Rule to manually remove duplicates present in samples.
rule remove_duplicates_manually:
    input:
        tsv='%s/fragments/fragmentomics/{sample}_{chr}_{pos}.tsv' % R_FOLDER,
    output:
        tsv='%s/fragments/remove_duplicates_manually/{sample}_{chr}_{pos}_{alt}_remove_duplicates.tsv' % R_FOLDER,
    benchmark:
        "%s/fragments/remove_duplicates_manually/{sample}_{chr}_{pos}_{alt}_remove_duplicates.tsv" % B_FOLDER
    log:
        "%s/fragments/remove_duplicates_manually/{sample}_{chr}_{pos}_{alt}_remove_duplicates.log" % L_FOLDER
    conda:
        "../envs/r.yaml"
    threads: 1
    resources:
        queue="shortq",
        mem_mb=500,
        time_min=20
    shell:
        """
        
        Rscript workflow/scripts/03.1_remove_duplicates_manually.R {input.tsv} {wildcards.sample} {wildcards.alt} {output.tsv} 
        
        """

