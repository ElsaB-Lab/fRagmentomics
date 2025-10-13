for file in *.R; do
  vim -c "%s/^  \(\S\)/    \1/g" -c "wq" "$file"
done
