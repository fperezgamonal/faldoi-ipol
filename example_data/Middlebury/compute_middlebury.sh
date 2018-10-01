#!/bin/bash
BIN_PATH=/home/fperez/Desktop/IPOL_article_figures/faldoi-ipol/scripts_python
SRC_PATH=/home/fperez/Desktop/IPOL_article_figures/faldoi-ipol/example_data/Middlebury
for f in *.txt; 
do  
    name=$(echo "$f" | cut -f 1 -d '.')
    echo "Processing folder $name..."
    echo "Filename: $f"
    cd $BIN_PATH
    time ./faldoi_deep.py $SRC_PATH/$f -nt 4 -vm 4 -split_img 0 -h_parts 3 -v_parts 2 -res_path ../Results/Middlebury/$name/tvl2csad_no_partitions/

done
