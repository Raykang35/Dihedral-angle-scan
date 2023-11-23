#!/bin/bash

# Check if a filename is provided as an argument
if [ $# -eq 0 ]; then
  echo "Usage: $0 <filename>"
  exit 1
fi

# Get the filename from the first argument
xyz="$1"

# Get only the name of the file
filename="${xyz%.*}"

# Take a look of xyz file
cat "$xyz"

# Create Gaussain input file
{ 
    echo "%mem=16GB"
    echo "%nprocshared=11"
    echo "#P b3lyp/6-31G* opt EmpiricalDispersion=GD3 gfinput pop=full"
    echo " "
    echo $filename
    echo " "
    echo "0 1"
} > "$filename.com"

# Grep the coordinates of xyz file and export to the com file
tail -n +3 "$xyz" | grep -E  "[a-zA-Z0-9]" >> "$filename.com"

# Add the blank line in the end of the com file
echo " " >> "$filename.com"

# Vim run 
vim "run"

# Exit
exit 0
