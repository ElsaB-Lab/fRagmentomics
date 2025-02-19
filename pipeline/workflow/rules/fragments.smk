# Extract reads covering small window around mutation of interest
# Then extract only high-quality reads
# Then select only first 11 columns
rule samtools_extract_reads:
    input:
        bam="%s/data/bam/{sample}.bam" % R_FOLDER
    output:
        sam_pos='%s/fragments/samtools_extract_reads/{sample}_{chr}_{pos}.sam' % R_FOLDER,
        sam_ext='%s/fragments/samtools_extract_reads/{sample}_{chr}_{pos}_extended.sam' % R_FOLDER
    benchmark:
        "%s/fragments/samtools_extract_reads/{sample}_{chr}_{pos}.tsv" % B_FOLDER
    log:
        "%s/fragments/samtools_extract_reads/{sample}_{chr}_{pos}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    params:
        pos_l=lambda w: int(w.pos) + config["params"]["fragments"]["neg_offset"],
        pos_r=lambda w: int(w.pos) + config["params"]["fragments"]["pos_offset"],
        flag_keep=config["params"]["picard_flag"]["flag_keep"],
        flag_remove=config["params"]["picard_flag"]["flag_remove"]
    threads: 1
    resources:
        queue="shortq",
        mem_mb=4000,
        time_min=90
    shell:
        """
        samtools view -b {input.bam} {wildcards.chr}:{wildcards.pos}-{wildcards.pos} | \
            samtools view -f {params.flag_keep} -F {params.flag_remove} - | \
            cut -f1-11 - > {output.sam_pos}

        samtools view -b {input.bam} {wildcards.chr}:{params.pos_l}-{params.pos_r} | \
            samtools view -f {params.flag_keep} -F {params.flag_remove} - | \
            cut -f1-11 - > {output.sam_ext}
       """


# Subset the extended SAM file to select only reads from fragments in the SAM file of reads
# covering the position of interest.
rule select_fragments_of_interest:
    input:
        sam_pos='%s/fragments/samtools_extract_reads/{sample}_{chr}_{pos}.sam' % R_FOLDER,
        sam_ext='%s/fragments/samtools_extract_reads/{sample}_{chr}_{pos}_extended.sam' % R_FOLDER
    output:
        sam='%s/fragments/select_fragments_of_interest/{sample}_{chr}_{pos}.sam' % R_FOLDER
    benchmark:
        "%s/fragments/select_fragments_of_interest/{sample}_{chr}_{pos}.tsv" % B_FOLDER
    log:
        "%s/fragments/select_fragments_of_interest/{sample}_{chr}_{pos}.log" % L_FOLDER
    threads: 1
    resources:
        queue="shortq",
        mem_mb=4000,
        time_min=90
    shell:
        """
        bash workflow/scripts/02.1_select_fragments_of_interest.sh \
            -a {input.sam_pos} \
            -b {input.sam_ext} \
            -o {output.sam} &> {log}
        """


# Core rule which extracts fragment statistics from the SAM of interest.
rule compute_fragmentomics:
    input:
        sam='%s/fragments/select_fragments_of_interest/{sample}_{chr}_{pos}.sam' % R_FOLDER,
        env='%s/setup/setup_r.done' % L_FOLDER
    output:
        tsv='%s/fragments/fragmentomics/{sample}_{chr}_{pos}.tsv' % R_FOLDER,
    benchmark:
        "%s/fragments/fragmentomics/{sample}_{chr}_{pos}.tsv" % B_FOLDER
    log:
        "%s/fragments/fragmentomics/{sample}_{chr}_{pos}.log" % L_FOLDER
    params:
        # Get ref, alt, mutation_type, del_info for the good combinaison of (sample, chr, pos)
        ref=lambda wc: get_ref(wc.sample, wc.chr, wc.pos),
        alt=lambda wc: get_alt(wc.sample, wc.chr, wc.pos),
        mutation_type=lambda wc: get_mutation_type(wc.sample, wc.chr, wc.pos),
        del_info=lambda wc : get_del_info(wc.sample, wc.chr, wc.pos)
    conda:
        "../envs/r.yaml"
    threads: 2
    resources:
        queue="shortq",
        mem_mb=30000,
        time_min=200
    shell:
        """
        Rscript workflow/scripts/02.2_compute_fragmentomics.R \
            --sam {input.sam} \
            --sample_id {wildcards.sample}
            --chr {wildcards.chr} \
            --pos {wildcards.pos} \
            --ref {params.ref} \
            --alt {params.alt} \
            --del_info {params.del_info} \
            --mutation_type {params.mutation_type} \
            --out {output.tsv} \
            --n_cores {threads} \
            --log {log}
        """



# Extract reads covering small window around mutation of interest
# Then extract only high-quality reads
# Then select only first 11 columns
#rule samtools_extract_reads:
#    input:
#        bam="%s/data/bam_without_dup/{sample}_without_dup.bam" % R_FOLDER
#    output:
#        sam_pos='%s/fragments/samtools_extract_reads/{sample}_{chr}_{pos}_without_dup.sam' % R_FOLDER,
#        sam_ext='%s/fragments/samtools_extract_reads/{sample}_{chr}_{pos}_without_dup_extended.sam' % R_FOLDER
#    benchmark:
#        "%s/fragments/samtools_extract_reads/{sample}_{chr}_{pos}_without_dup.tsv" % B_FOLDER
#    log:
#        "%s/fragments/samtools_extract_reads/{sample}_{chr}_{pos}_without_dup.log" % L_FOLDER
#    conda:
#        "../envs/main.yaml"
#    params:
#        pos_l=lambda w: int(w.pos) + config["params"]["fragments"]["neg_offset"],
#        pos_r=lambda w: int(w.pos) + config["params"]["fragments"]["pos_offset"]
#    threads: 1
#    resources:
#        queue="shortq",
#        mem_mb=4000,
#        time_min=90
#    shell:
#        """
#        samtools view -b {input.bam} {wildcards.chr}:{wildcards.pos}-{wildcards.pos} | \
#            samtools view -f 0x03 -F 0x900 - | \
#            cut -f1-11 - > {output.sam_pos}
#        samtools view -b {input.bam} {wildcards.chr}:{params.pos_l}-{params.pos_r} | \
#            samtools view -f 0x03 -F 0x900 - | \
#            cut -f1-11 - > {output.sam_ext}
#        """

# Subset the extended SAM file to select only reads from fragments in the SAM file of reads
# covering the position of interest.
#rule select_fragments_of_interest:
#    input:
#        sam_pos='%s/fragments/samtools_extract_reads/{sample}_{chr}_{pos}_without_dup.sam' % R_FOLDER,
#        sam_ext='%s/fragments/samtools_extract_reads/{sample}_{chr}_{pos}_without_dup_extended.sam' % R_FOLDER
#    output:
#        sam='%s/fragments/select_fragments_of_interest/{sample}_{chr}_{pos}_without_dup.sam' % R_FOLDER
#    benchmark:
#        "%s/fragments/select_fragments_of_interest/{sample}_{chr}_{pos}_without_dup.tsv" % B_FOLDER
#    log:
#        "%s/fragments/select_fragments_of_interest/{sample}_{chr}_{pos}_without_dup.log" % L_FOLDER
#    threads: 1
#    resources:
#        queue="shortq",
#        mem_mb=4000,
#        time_min=90
#    shell:
#        """
#        bash workflow/scripts/02.1_select_fragments_of_interest.sh \
#            -a {input.sam_pos} \
#            -b {input.sam_ext} \
#            -o {output.sam} &> {log}
#        """


# Core rule which extracts fragment statistics from the SAM of interest.
#rule compute_fragmentomics:
#    input:
#        sam='%s/fragments/select_fragments_of_interest/{sample}_{chr}_{pos}_without_dup.sam' % R_FOLDER,
#        env='%s/setup/setup_r.done' % L_FOLDER
#    output:
#        tsv='%s/fragments/fragmentomics/{sample}_{chr}_{pos}_without_dup.tsv' % R_FOLDER,
#    benchmark:
#        "%s/fragments/fragmentomics/{sample}_{chr}_{pos}_without_dup.tsv" % B_FOLDER
#    log:
#        "%s/fragments/fragmentomics/{sample}_{chr}_{pos}_without_dup.log" % L_FOLDER
#    conda:
#        "../envs/r.yaml"
#    threads: 20
#    resources:
#        queue="shortq",
#        mem_mb=30000,
#        time_min=200
#    shell:
#        """
#        Rscript workflow/scripts/02.2_compute_fragmentomics.R \
#            --sam {input.sam} \
#            --chr {wildcards.chr} \
#            --pos {wildcards.pos} \
#            --out {output.tsv} \
#            --n_cores {threads} \
#            --log {log}
#        """


