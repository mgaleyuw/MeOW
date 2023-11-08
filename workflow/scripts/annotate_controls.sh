#!/bin/bash

#merge output of meOW and join to useful annotations, optionally. Take a directory or a list of tsv files to merge.

source /home/mgaley/methScripts/defaultConfig.config

LC_ALL=C

fileList=( ${snakemake_input[@]} )
outputName=${snakemake_output[0]}

if [[ -f ${snakemake_params[startingFile]} ]]

if [ "${#fileList[@]}" -lt 2 ]
then
	echo "need a minimum of 2 files to run join"
	exit 1
fi

if [[ -f ${snakemake_params[startingFile]} ]]
then
    fileList=( ${snakemake_params[startingFile]} ${fileList[@]} )
fi

firstFile=${fileList[0]}


for ifile in ${fileList[@]:1}
do
	secondFile=$ifile
	numFields=$( awk '{print NF;exit}' $firstFile)
	mergeString1=$(for i in $(seq 2 $(echo "$numFields" | bc)); do echo "1.$i"; done | tr '\n' ' ')
	echo "$numFields"
	numFields=$( awk '{print NF;exit}' $secondFile)
	mergeString2=$(for i in $(seq 2 $(echo "$numFields" | bc)); do echo "2.$i"; done | tr '\n' ' ')
	
	fullMergeString=$(echo "0 $mergeString1$mergeString2" | xargs)
	echo "'"$fullMergeString"'"
	LC_ALL=C
	join <(head -n 1 $firstFile) <(head -n 1 $secondFile) -t $'\t' > mergedFile.tsv
	LC_ALL=C
	join <(tail -n+2 $firstFile | sort -k 1,1) <(tail -n+2 $secondFile | sort -k 1,1) -t $'\t' -e 'NA' -a1 -a2 -o "$fullMergeString" >> mergedFile.tsv
	cp mergedFile.tsv newFirstFile.tsv
	firstFile=newFirstFile.tsv
done

rm newFirstFile.tsv

#mergedFile.tsv is the input for annotation

INPUT=mergedFile.tsv
# the input file has two columns but would work with more. expects a header, which should be saved. 

echo -e "chr.pos\tcpg" > intermediate.names.tsv
cat $INPUT | sed 's/\./\t/' |tail -n+2 | sort -k1,1 -k2,2n | cut -f 1,2 | $RUSTANNOTATE $NAMEDCPGS h | sort -k1,1 >> intermediate.names.tsv

echo -e "$(head -n 1 $INPUT)\tcpg" > intermediate.cpgs.named.all.tsv
join <(tail -n+2 $INPUT | sort -k1,1) <(tail -n+2 intermediate.names.tsv | sort -k1,1) -t $'\t' >> intermediate.cpgs.named.all.tsv
# produces a two column file 

#get gene names
cat $INPUT | sed 's/\./\t/' |tail -n+2 | sort -k1,1 -k2,2n | cut -f 1,2 | $RUSTANNOTATE $NAMEDGENES h | sort -k1,1 > intermediate.gene.names.tsv

echo -e "$(head -n 1 intermediate.cpgs.named.all.tsv)\tgene" > cpgs.named.all.tsv
join <(tail -n+2 intermediate.cpgs.named.all.tsv | sort -k1,1) <(tail -n+2 intermediate.gene.names.tsv | sort -k1,1) -t $'\t' >> cpgs.named.all.tsv

#finally fix the indexing
paste <(cat cpgs.named.all.tsv | cut -f 1 | tr '.' '\t') <(cut -f2- cpgs.named.all.tsv) > $outputName

rm intermediate*tsv
rm cpgs.named.all.tsv
