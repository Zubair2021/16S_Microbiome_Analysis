#!/bin/bash

function get_all_files {
if [ ! -d $1 ]; then
echo "Error: $1 is not a valid directory."
return
fi

parent_dir=$(dirname $1)
new_dir="$parent_dir/fastq_files"
mkdir $new_dir

find $1 -type f -exec cp {} $new_dir \;
}

get_all_files $1
