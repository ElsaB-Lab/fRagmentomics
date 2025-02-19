#!/bin/bash

while getopts ":o:" opt; do
    case $opt in
	o) output_file="$OPTARG"
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

# Get output file basename
output_basename=$(basename $output_file)

# Extract Job ID and Node from the output file basename
job_id=$(echo "$output_basename" | grep -oP '(?<=_)\d+(?=-n)')
node=$(echo "$output_basename" | grep -oP '(?<=n)\d+')

# Extract the input .sam file path
sam_file=$(grep -oP 'input: \K[^ ]+\.sam' $output_file)

# Capture file size and number of lines
sam_file_size=$(stat -c "%s" "$sam_file")
sam_num_lines=$(wc -l < "$sam_file")

# Extract relevant information from the job output
sample_id=$(grep -oP 'sample=\K[^,]+' $output_file)
chromosome=$(grep -oP 'chr=\K[^,]+' $output_file)
position=$(grep -oP 'pos=\K[^,]+' $output_file)
rule_name=$(grep -oP 'rule \K\w+' $output_file)

# Extract CPU Requested, Memory Requested, and Time Requested from the job output
cpu_requested=$(grep -oP '^\s*threads: \K\d+' $output_file)
memory_requested=$(grep -oP '^\s*resources:.*\bmem_mb=\K\d+' $output_file)
time_requested=$(grep -oP '^\s*resources:.*\btime_min=\K\d+' $output_file)

# Get seff output
seff_output=$(seff $job_id)

# Extract CPU Efficiency, Memory Efficiency, and Wall Clock Time from the seff output
cpu_efficiency=$(echo "$seff_output" | awk '/CPU Efficiency:/ {print $3}')
memory_efficiency=$(echo "$seff_output" | awk '/Memory Efficiency:/ {print $3}')
clock_=$(echo "$seff_output" | awk '/Memory Efficiency:/ {print $3}')
wall_clock_time=$(echo "$seff_output" | grep -oP 'Job Wall-clock time: \K[0-9][0-9]:[0-9][0-9]:[0-9][0-9]')
wall_clock_time=$(echo "$wall_clock_time" | awk -F: '{print $1 * 60 + $2}')

# Extract the benchmark .tsv file path  
tsv_file=$(grep -oP 'benchmark: \K[^ ]+\.tsv' $output_file) 

# Capture file s and max_vms 
running_time_sec=$(awk 'FNR == 2 {print $1}' "$tsv_file") 
max_vms=$(awk 'FNR == 2 {print $4}' "$tsv_file") 

# Print the output in tab-separated format
printf "Sample_Id\tChromosome\tPosition\tRule_Name\tJob_Id\tNode\tSam_File_Size\tSam_Num_Lines\tCPU_Requested\tCPU_Efficiency\tMemory_Requested\tMemory_Efficiency\tTime_Requested\tWall_Clock_Time_Min\tRunning_Time_Sec\tMax_Vms\n"
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$sample_id" "$chromosome" "$position" "$rule_name" "$job_id" "$node" "$sam_file_size" "$sam_num_lines" "$cpu_requested" "$cpu_efficiency" "$memory_requested" "$memory_efficiency" "$time_requested" "$wall_clock_time" "$running_time_sec" "$max_vms"
