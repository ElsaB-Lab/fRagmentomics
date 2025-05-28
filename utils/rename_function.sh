#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <old_function_name> <new_function_name>"
    exit 1
fi

OLD_NAME=$1
NEW_NAME=$2

# Rename function in .R files in the R directory
find R -type f -name "*.R" -exec sed -i '' "s/$OLD_NAME/$NEW_NAME/g" {} \;

# Rename function in test files in the tests/testthat directory
find tests/testthat -type f -name "*.R" -exec sed -i '' "s/$OLD_NAME/$NEW_NAME/g" {} \;

# Rename filenames in the R directory
find R -type f -name "*$OLD_NAME*.R" | while read -r file; do
    new_file=$(echo "$file" | sed "s/$OLD_NAME/$NEW_NAME/")
    mv "$file" "$new_file"
done

# Rename filenames in the tests/testthat directory
find tests/testthat -type f -name "*$OLD_NAME*.R" | while read -r file; do
    new_file=$(echo "$file" | sed "s/$OLD_NAME/$NEW_NAME/")
    mv "$file" "$new_file"
done

echo "Function and file renaming complete."
