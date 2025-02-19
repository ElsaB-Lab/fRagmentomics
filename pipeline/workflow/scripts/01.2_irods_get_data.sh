#!/bin/bash

while getopts ":r:s:t:a:b:c:d:w:x:y:z:" opt; do
    case $opt in
	r) dir_bam="$OPTARG"
	    ;;
	s) dir_pdf="$OPTARG"
	    ;;
	t) dir_xml="$OPTARG"
	    ;;
	a) irods_bam="$OPTARG"
	    ;;
	b) irods_bai="$OPTARG"
	    ;;
	c) irods_pdf="$OPTARG"
	    ;;
	d) irods_xml="$OPTARG"
	    ;;
	w) local_bam="$OPTARG"
	    ;;
	x) local_bai="$OPTARG"
	    ;;
	y) local_pdf="$OPTARG"
	    ;;
	z) local_xml="$OPTARG"
	    ;;
	\?) echo "Invalid option -$OPTARG" >&2
	    exit 1
	    ;;
    esac

    case $OPTARG in
	-*) echo "Option $opt needs a valid argument"
	    exit 1
	    ;;
    esac
done

basename_bam="$(basename $irods_bam)"
basename_bai="$(basename $irods_bai)"
basename_pdf="$(basename $irods_pdf)"
basename_xml="$(basename $irods_xml)"

# get BAM file and rename if necessary
mkdir -p $dir_bam
iget -vK ${irods_bam} ${dir_bam}

if [[ "${dir_bam}/${basename_bam}" != "${local_bam}" ]]; then
    mv ${dir_bam}/${basename_bam} ${local_bam}
    echo "moved ${dir_bam}/${basename_bam} to ${local_bam}!"
fi

# get BAI file and rename if necessary
mkdir -p $dir_bam
iget -vK ${irods_bai} ${dir_bam}

if [[ "${dir_bam}/${basename_bai}" != "${local_bai}" ]]; then
    mv ${dir_bam}/${basename_bai} ${local_bai}
    echo "moved ${dir_bam}/${basename_bai} to ${local_bai}!"
fi

# get PDF file and rename if necessary
mkdir -p $dir_pdf
iget -vK ${irods_pdf} ${dir_pdf}

if [[ "${dir_pdf}/${basename_pdf}" != "${local_pdf}" ]]; then
    mv ${dir_pdf}/${basename_pdf} ${local_pdf}
    echo "moved ${dir_pdf}/${basename_pdf} to ${local_pdf}!"
fi

# get XML file and rename if necessary
mkdir -p $dir_xml
iget -vK ${irods_xml} ${dir_xml}

if [[ "${dir_xml}/${basename_xml}" != "${local_xml}" ]]; then
    mv ${dir_xml}/${basename_xml} ${local_xml}
    echo "moved ${dir_xml}/${basename_xml} to ${local_xml}!"
fi
