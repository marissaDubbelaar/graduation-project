#!/bin/bash
# Before calling this file you must go to the directory with all of the .counts files
# that you want to merge.
#
# The files that end with .count in this directory will be merged together into one
# file. This script is made because some of the R packages for the RNA-sequencing analysis
# prefer that the data is in one file.

############################################################
#							Code						   #
############################################################
# This pattern gets all of the filenames with that end with .count
FILES=$(ls -t -v *.txt | tr '\n' ' ');
# The code below makes sure that the header for the data set is made
fileNames="probes"
for i in $FILES;
do

	fileNames="$fileNames\t$(echo $i | awk -F '[.]' '{print $1}')"
done
# The filenames will be put into the file first.
# Creating the header of the mergeCounts.txt
echo $fileNames > mergedCountstmp.txt
# The code below writes the output with the use of the AWK language.
# The NF contains the number of field within the file.
# If the NF is more than one, the files will be put next to each other (tab separated).
# The files found in the directory will be compared one by one, the results are seved into
# mergedCounts.txt
awk 'NF > 1{ a[$1] = a[$1]"\t"$2} END {for( i in a ) print i a[i]}' $FILES >> mergedCountstmp.txt
# Remove the rows that start with '__'
awk '!/^__/' mergedCountstmp.txt > mergedCounts.txt
rm mergedCountstmp.txt
