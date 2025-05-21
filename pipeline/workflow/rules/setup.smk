rule setup_r:
    input:
        f"../fRagmentomics_{config['tools_versions']['fragmentomics']}.tar.gz"
    output:
        touch("%s/setup/setup_r.done" % L_FOLDER)
    conda:
        "../envs/r.yaml"
    log:
        "%s/setup/setup_r.log" % L_FOLDER
    benchmark:
        "%s/setup/setup_r.tsv" % B_FOLDER
    resources:
        queue="shortq",
        mem_mb=16000,
        time_min=60
    shell:
        """
        R CMD INSTALL {input}
        """
