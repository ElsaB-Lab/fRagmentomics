if not config["general"]["automatic_irods_check"]["disable"]:
    rule irods_get_content:
        log:
            "%s/irods/irods_get_content.log" % L_FOLDER
        benchmark:
            "%s/irods/irods_get_content.log" % B_FOLDER
        params:
            project_names = config["general"]["automatic_irods_check"]["project_names"]
        output:
            irods = "%s/irods/irods_get_content.tsv" % R_FOLDER
        resources:
            queue = "shortq"
        threads: 1
        shell:
            """
            module unload anaconda3
            pns=( {params.project_names} )
            for pn in "${{pns[@]}}"; do
                imeta qu -C projectName like "%${{pn}}%" >> {output.irods} 2> {log}
            done
            """

    rule irods_check_status:
        log:
            "%s/irods/irods_check_status.log" % L_FOLDER
        benchmark:
            "%s/irods/irods_check_status.log" % B_FOLDER
        conda:
            "../envs/python.yaml"
        input:
            table = config["general"]["samples"],
            irods = "%s/irods/irods_get_content.tsv" % R_FOLDER
        output:
            status = "%s/irods/irods_check_status.pass" % R_FOLDER
        resources:
            queue = "shortq"
        threads: 1
        shell:
            """
            python workflow/scripts/01.1_irods_check_status.py \
                --table {input.table} \
                --irods {input.irods} \
                --output {output.status} &> {log}
            """
else:
    open("%s/irods/irods_check_status.pass" % R_FOLDER, "w").close()


rule irods_get_data:
    log:
        "%s/irods/irods_get_data_{sample}.log" % L_FOLDER
    benchmark:
        "%s/irods/irods_get_data_{sample}.log" % B_FOLDER
    input:
        table = config["general"]["samples"],
        check = "%s/irods/irods_check_status.pass" % R_FOLDER
    params:
        dir_bam = "%s/data/bam" % R_FOLDER,
        dir_pdf = "%s/data/pdf" % R_FOLDER,
        dir_xml = "%s/data/xml" % R_FOLDER,
        irods_bam = lambda w: get_data_irods(w)["BAM"],
        irods_bai = lambda w: get_data_irods(w)["BAI"],
        irods_pdf = lambda w: get_data_irods(w)["PDF"],
        irods_xml = lambda w: get_data_irods(w)["XML"]
    output:
        bam = "%s/data/bam/{sample}.bam" % R_FOLDER,
        bai = "%s/data/bam/{sample}.bam.bai" % R_FOLDER,
        pdf = "%s/data/pdf/{sample}.pdf" % R_FOLDER,
        xml = "%s/data/xml/{sample}.xml" % R_FOLDER
    resources:
        local_load = 1
    threads: 1
    shell:
        """
        module unload anaconda3
        bash workflow/scripts/01.2_irods_get_data.sh \
            -r {params.dir_bam} \
            -s {params.dir_pdf} \
            -t {params.dir_xml} \
            -a {params.irods_bam} \
            -b {params.irods_bai} \
            -c {params.irods_pdf} \
            -d {params.irods_xml} \
            -w {output.bam} \
            -x {output.bai} \
            -y {output.pdf} \
            -z {output.xml} &> {log}
        """
