#!/bin/bash

### Create scan rigid input files file ###

# Check if a filename is provided as an argument
if [ $# -eq 1 ]; then
  echo "Usage: $0 <filename>, <Diangle#>"
  exit 1
fi

# Get the filename from the first argument
xyz="$1"
diangle="$2"

# Get only the name of the file
filename="${xyz%.*}"
Diangle="${diangle%.*}"

# Create Gaussain input file
{
    echo "%mem=16GB"
    echo "%nprocshared=11"
    echo "#P b3lyp/6-31G* EmpiricalDispersion=GD3 gfinput pop=full"
    echo " "
    echo "$filename"
    echo " "
    echo "0 1"
} > "$diangle.com"

# Grep the coordinates of xyz file and export to the com file
tail -n +3 "$xyz" | grep -E  "[a-zA-Z0-9]" >> "$diangle.com"

# Add the blank line in the end of the com file
echo " " >> "$diangle.com"

# Fill the filename
sed "s/dxxx/${Diangle}/g" run >> run_${Diangle}; rm run

# Exit
exit 0
