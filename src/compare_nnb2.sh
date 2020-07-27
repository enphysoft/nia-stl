#!/bin/bash
#+TITLE: compare.sh
#+File: /home/albertsk/Documents/research/publications/KRISO-Cball-STL-nnbor/codes/nia-stl/src/compare.sh
#+Date: Fri Jul 24 15:07:21 HST 2020
#+Author: Albert S. Kim, albertsk@hawaii.edu
#
#
. ~/.bashrc 
# . /opt/openfoam5/etc/bashrc
file1=$1
file2=$2
# STL_INPUT_nia1.stl		# 
# cube0-from-gmsh_nnbs.stl 

# grep endfacet $file1 | awk -F, 'BEGIN { OFS = ","} {print $1,$2,$3,$4,$5,$6," "}' | sed "s/endfacet//g" | sed "s/ //g" > temp1.dat
# grep endfacet $file2 | awk -F, 'BEGIN { OFS = ","} {print $1,$2,$3,$4, "   "}' | sed "s/endfacet//g" | sed "s/ //g" > temp2.dat
grep endfacet $file1 | sed "s/endfacet//g" | sed "s/ //g" > temp1.dat
grep endfacet $file2 | sed "s/endfacet//g" | sed "s/ //g" > temp2.dat

paste temp1.dat temp2.dat
diff  temp1.dat temp2.dat

# diff temp1.dat temp2.dat
# rm -f temp1.dat temp2.dat

numlines=$(cat temp1.dat | wc -l)
echo 
echo " Comparing neighbors of two output fils of" $numlines "lines."
for i in `seq 1 $numlines`
do
    echo "=== Compare Facet" $i "of <"$file1"> and <"$file2">"  
    sed -n "${i}p" temp1.dat
    sed -n "${i}p" temp2.dat
    echo 
done



exit 0



