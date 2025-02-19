# Run FASTQC on BAM files
# See https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
rule fastqc:
    input:
        bam="%s/data/bam/{sample}.bam" % R_FOLDER,
        bai="%s/data/bam/{sample}.bam.bai" % R_FOLDER
    output:
        html='%s/qc/fastqc/{sample}_fastqc.html' % R_FOLDER,
        zip='%s/qc/fastqc/{sample}_fastqc.zip' % R_FOLDER
    benchmark:
        "%s/qc/fastqc/{sample}.tsv" % B_FOLDER
    log:
        "%s/qc/fastqc/{sample}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    params:
        out="%s/qc/fastqc" % R_FOLDER
    threads: 1
    resources:
        queue="mediumq",
        mem_mb=5000,
        time_min=600
    shell:
        'fastqc -o {params.out} {input.bam} 2> {log}'


# Quality-control of mapping with samtools stats
# See http://www.htslib.org/doc/samtools-stats.html
rule samtools_stats:
    input:
        bam="%s/data/bam/{sample}.bam" % R_FOLDER,
        bai="%s/data/bam/{sample}.bam.bai" % R_FOLDER
    output:
        "%s/qc/samtools_stats/{sample}_stats.tsv" % R_FOLDER
    benchmark:
        "%s/qc/samtools_stats/{sample}.tsv" % B_FOLDER
    log:
        "%s/qc/samtools_stats/{sample}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    threads: 1
    resources:
        queue="shortq",
        mem_mb=3000,
        time_min=60
    shell:
        """samtools stats {input.bam} > {output} 2> {log}"""


# Quality-control of mapping with samtools flagstat
# See http://www.htslib.org/doc/samtools-flagstat.html
rule samtools_flagstat:
    input:
        bam="%s/data/bam/{sample}.bam" % R_FOLDER,
        bai="%s/data/bam/{sample}.bam.bai" % R_FOLDER
    output:
        "%s/qc/samtools_flagstat/{sample}_flagstat.tsv" % R_FOLDER
    benchmark:
        "%s/qc/samtools_flagstat/{sample}.tsv" % B_FOLDER
    log:
        "%s/qc/samtools_flagstat/{sample}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    threads: 1
    resources:
        queue="shortq",
        mem_mb=3000,
        time_min=60
    shell:
        """samtools flagstat {input.bam} > {output} 2> {log}"""



# Quality-control of mapping with gatk CollectMultipleMetrics (picard)
# See https://gatk.broadinstitute.org/hc/en-us/articles/360037594031-CollectMultipleMetrics-Picard
rule collect_multiple_metrics:
    input:
        bam="%s/data/bam/{sample}.bam" % R_FOLDER,
        bai="%s/data/bam/{sample}.bam.bai" % R_FOLDER
    output:
        directory("%s/qc/collect_multiple_metrics/collect_multiple_metrics_{sample}" % R_FOLDER)
    benchmark:
        "%s/qc/collect_multiple_metrics/{sample}.tsv" % B_FOLDER
    log:
        "%s/qc/collect_multiple_metrics/{sample}.log" % L_FOLDER
    conda:
         "../envs/main.yaml"
    params:
        java="'-Xmx40g -Djava.io.tmpdir=%s'" % config["general"]["tmpdir"],
    threads: 1
    resources:
        queue="shortq",
        mem_mb=30000,
        time_min=90
    shell:
        """gatk --java-options {params.java} CollectMultipleMetrics \
            --I {input.bam} \
            --PROGRAM CollectAlignmentSummaryMetrics \
            --PROGRAM CollectInsertSizeMetrics \
            --PROGRAM CollectSequencingArtifactMetrics \
            --O {output} 2> {log}"""




# Remove duplicated reads from the bam file with MarkDuplicates (picard)
# See https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard
rule mark_duplicates:
    input:
        "%s/data/bam/{sample}.bam" % R_FOLDER
    output:
        bam="%s/qc/mark_duplicates/{sample}_without_dup.bam" % R_FOLDER,
        metrics="%s/qc/mark_duplicates/{sample}.metrics.txt" % R_FOLDER
    benchmark:
        "%s/qc/mark_duplicates/{sample}.tsv" % B_FOLDER
    log:
        "%s/qc/mark_duplicates/{sample}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    params:
        java="'-Xmx40g -Djava.io.tmpdir=%s'" % config["general"]["tmpdir"],
    threads: 1
    resources:
        queue="shortq",
        mem_mb=50000,
        time_min=240
    shell:
        """gatk --java-options {params.java} MarkDuplicates \
            --I {input} \
            --REMOVE_DUPLICATES {params.remove_duplicates} \
            --O {output.bam} \
            --M {output.metrics} 2> {log}"""



rule multiqc:
    input:
        expand("%s/qc/fastqc/{sample}_fastqc.zip" % R_FOLDER, sample=samples),
        expand("%s/qc/samtools_flagstat/{sample}_flagstat.tsv" % R_FOLDER, sample=samples),
        expand("%s/qc/samtools_stats/{sample}_stats.tsv" % R_FOLDER, sample=samples),
        expand("%s/qc/mark_duplicates/{sample}.metrics.txt" % R_FOLDER, sample=samples),
        expand("%s/qc/collect_multiple_metrics/collect_multiple_metrics_{sample}.alignment_summary_metrics" % R_FOLDER, sample=samples),
        expand("%s/qc/collect_multiple_metrics/collect_multiple_metrics_{sample}.insert_size_metrics" % R_FOLDER, sample=samples),
        expand("%s/qc/collect_multiple_metrics/collect_multiple_metrics_{sample}.base_distribution_by_cycle_metrics" % R_FOLDER, sample=samples),
        expand("%s/qc/collect_multiple_metrics/collect_multiple_metrics_{sample}.quality_by_cycle_metrics" % R_FOLDER, sample=samples),
        expand("%s/qc/collect_multiple_metrics/collect_multiple_metrics_{sample}.quality_distribution_metrics" % R_FOLDER, sample=samples)
    output:
        '%s/qc/multiqc/multiqc_report.html' % R_FOLDER,
        directory('%s/qc/multiqc/multiqc_data' % R_FOLDER)
    benchmark:
        "%s/qc/multiqc/multiqc.tsv" % B_FOLDER
    log:
        "%s/qc/multiqc/multiqc.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    params:
        inpdir_1="%s/qc/fastqc" % R_FOLDER,
        inpdir_2="%s/qc/samtools_flagstat" % R_FOLDER,
        inpdir_3="%s/qc/samtools_stats" % R_FOLDER,
        inpdir_4="%s/qc/mark_duplicates" % R_FOLDER,
        inpdir_5="%s/qc/collect_multiple_metrics" % R_FOLDER,
        outdir="%s/qc/multiqc" % R_FOLDER
    threads: 1
    resources:
        queue="shortq",
        mem_mb=16000,
        time_min=30
    shell:
        'multiqc {params.inpdir_1} {params.inpdir_2} {params.inpdir_3} {params.inpdir_4} {params.inpdir_5} --outdir {params.outdir} &> {log}'



# # Quality-control of mapping with gatk CollectHsMetrics (picard)
# # See https://gatk.broadinstitute.org/hc/en-us/articles/360036856051-CollectHsMetrics-Picard-
# rule collect_hs_metrics:
#     input:
#         bam="%s/data/bam/{sample}.bam" % R_FOLDER,
#         bai="%s/data/bam/{sample}.bam.bai" % R_FOLDER,
#         intervals=# path to target file
#     output:
#         "%s/qc/collect_hs_metrics/{sample}_hs_metrics.tsv" % R_FOLDER
#     benchmark:
#         "%s/qc/collect_hs_metrics/{sample}.tsv" % B_FOLDER
#     log:
#         "%s/qc/collect_hs_metrics/{sample}.log" % L_FOLDER
#     conda:
#         "../envs/main.yaml"
#     params:
#         java="'-Xmx40g -Djava.io.tmpdir=%s'" % config["general"]["tmpdir"],
#     threads: 1
#     resources:
#         queue="shortq",
#         mem_mb=30000,
#         time_min=90
#     shell:
#         """gatk --java-options {params.java} CollectHsMetrics \
#             --TI {input.intervals} \
#             --BI {input.intervals} \
#             --I {input.bam} \
#             --O {output} 2> {log}"""


# # Quality-control of mapping with mosdepth
# # See https://github.com/brentp/mosdepth
# rule mosdepth:
#     input:
#         bam="%s/data/bam/{sample}.bam" % R_FOLDER,
#         bai="%s/data/bam/{sample}.bam.bai" % R_FOLDER,
#     output:
#         "%s/qc/mosdepth/{sample}.mosdepth.global.dist.txt" % R_FOLDER
#     benchmark:
#         "%s/qc/mosdepth/{sample}.tsv" % B_FOLDER
#     log:
#         "%s/qc/mosdepth/{sample}.log" % L_FOLDER
#     conda:
#         "../envs/main.yaml"
#     params:
#         prefix="%s/qc/mosdepth/{sample}" % R_FOLDER
#     threads: 4
#     resources:
#         queue="shortq",
#         mem_mb=3000,
#         time_min=60
#     shell:
#         """mosdepth --no-per-base \
#             --threads{threads} \
#             --fast-mode \
#             --by {input.bed} \
#             {params.prefix} {input.bam} 2> {log}"""


# # Check BAM coverage
# rule bam_coverage:
#     input:
#         bam="%s/mapping/{sample}.nodup.recal.bam" % R_FOLDER,
#         bai="%s/mapping/{sample}.nodup.recal.bam.bai" % R_FOLDER
#     output:
#         "%s/qc/bam_coverage/{region}/{sample}.bw" % R_FOLDER
#     benchmark:
#         "%s/qc/bam_coverage/{sample}_{region}.tsv" % B_FOLDER
#     log:
#         "%s/qc/bam_coverage/{sample}_{region}.log" % L_FOLDER
#     conda:
#         "../envs/deeptools.yaml"
#     threads: 1
#     resources:
#         queue="shortq",
#         mem_mb=4000,
#         time_min=20
#     shell:
#         """
#         bamCoverage -b {input.bam} \
#             -r {wildcards.region} \
#             -bs 10 \
#             -o {output} &> {log}
#         """
