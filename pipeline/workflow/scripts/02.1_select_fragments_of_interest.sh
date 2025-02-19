#!/bin/bash

while getopts ":a:b:o:" opt; do
    case $opt in
	a) sam_pos="$OPTARG"
	    ;;
	b) sam_ext="$OPTARG"
	    ;;
	o) sam_out="$OPTARG"
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

# Extract the fragment names of interest from the first column ${sam_pos}
fragments_of_interest=$(cut -f 1 ${sam_pos} | sort | uniq)

# Use grep to filter the lines in ${sam_ext} that match any fragment of interest
grep -Fwf <(echo "$fragments_of_interest") ${sam_ext} > ${sam_out}
