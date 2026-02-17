#!/bin/bash

# Array to hold function names
declare -a functions

# Find all .R files in the current directory and subdirectories
while IFS= read -r file; do
  # Use grep to find function definitions and extract function names
  while IFS= read -r func; do
    functions+=("$func")
    echo "Found function: $func" # Debugging line to print each function name
  done < <(grep -E '^[[:space:]]*[a-zA-Z0-9._]+[[:space:]]+<-[[:space:]]+function[[:space:]]*\(' "$file" |
    sed -E 's/^[[:space:]]*([a-zA-Z0-9._]+)[[:space:]]+<-[[:space:]]+function[[:space:]]*\(.*$/\1/')
done < <(find . -type f -name "*.R")

# Sort the functions alphabetically
IFS=$'\n' sorted_functions=($(sort <<<"${functions[*]}"))
unset IFS

# Print sorted functions
echo "Sorted function names:"
for func in "${sorted_functions[@]}"; do
  echo "$func"
done

# Check for function names that are strictly included in others
echo -e "\nChecking for function names that are strictly included in others..."
for ((i = 0; i < ${#sorted_functions[@]}; i++)); do
  for ((j = i + 1; j < ${#sorted_functions[@]}; j++)); do
    if [[ "${sorted_functions[i]}" == *"${sorted_functions[j]}"* && "${sorted_functions[i]}" != "${sorted_functions[j]}" ]]; then
      echo "Warning: '${sorted_functions[j]}' is strictly included in '${sorted_functions[i]}"
    elif [[ "${sorted_functions[j]}" == *"${sorted_functions[i]}"* && "${sorted_functions[i]}" != "${sorted_functions[j]}" ]]; then
      echo "Warning: '${sorted_functions[i]}' is strictly included in '${sorted_functions[j]}"
    fi
  done
done
