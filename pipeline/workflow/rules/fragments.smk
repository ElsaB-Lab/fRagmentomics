# Select mutations for the sample
rule select_mutations:
    input:
        config["general"]["mutations"]
    output:
        '%s/data/mut/{sample}.tsv' % R_FOLDER
    benchmark:
        "%s/select_mutations/{sample}.tsv" % B_FOLDER
    log:
        "%s/select_mutations/{sample}.log" % L_FOLDER
    resources:
        queue="shortq",
        mem_mb=512,
        time_min=10
    shell:
        """
        awk -F'\\t' 'NR == 1 || $1 == "{wildcards.sample}"' {input} > {output} 2> {log}
        """


# Core rule which extracts fragment statistics from the SAM of interest.
rule fragmentomics:
    input:
        bam="%s/data/bam/{sample}.bam" % R_FOLDER,
        mut="%s/data/mut/{sample}.tsv" % R_FOLDER,
        fasta=config["ref"]["fasta"],
        env="%s/setup/setup_r.done" % L_FOLDER
    output:
        '%s/fragmentomics/{sample}.tsv' % R_FOLDER
    benchmark:
        "%s/fragmentomics/{sample}.tsv" % B_FOLDER
    log:
        "%s/fragmentomics/{sample}.log" % L_FOLDER
    params:
        neg_offset=config["params"]["fragments"]["neg_offset"],
        pos_offset=config["params"]["fragments"]["pos_offset"],
        flag_keep=config["params"]["fragments"]["flag_keep"],
        flag_remove=config["params"]["fragments"]["flag_remove"]
    threads: 8
    resources:
        queue="shortq",
        mem_mb=8000,
        time_min=120
    shell:
        """
        eval "$(conda shell.bash hook)"
        conda activate r-fragmentomics

        Rscript -e 'fRagmentomics::process_fragmentomics(
            mut="{input.mut}",
            bam="{input.bam}",
            fasta="{input.fasta}",
            sample_id="{wildcards.sample}",
            neg_offset_mate_search={params.neg_offset},
            pos_offset_mate_search={params.pos_offset},
            one_based=TRUE,
            flag_keep={params.flag_keep},
            flag_remove={params.flag_remove},
            report_tlen=TRUE,
            report_softclip=TRUE,
            report_5p_3p_bases_fragment=5,
            output_file="{output}",
            tmp_folder="/mnt/beegfs02/scratch/y_pradat/tmp",
            n_cores={threads}
        )' >> {log}
        """


