#!/bin/bash
trap 'trap " " SIGTERM; kill 0; wait' SIGINT SIGTERM EXIT

LC_ALL=C
CpGFile=${snakemake_config[DEFAULTREGION]}

SAMPLE=${snakemake_wildcards[SAMPLE]}
OUTPUTNAME=${snakemake_output[0]}
INPUTBAM=${snakemake_input[0]}
PARALLEL=${snakemake_params[parallel]}

REGIONNAMES=($(cut -f1 $CpGFile))
CHROMS=($(cut -f2 $CpGFile))
REGIONSTARTS=($(cut -f3 $CpGFile))
REGIONSTOPS=($(cut -f4 $CpGFile))


extractChromosome () {
  local REGIONNAME=$1
  local CHROM=$2
  local SAMPLE=$3
  local BAMFILE=$4

  if [[ -f "$BAMFILE.bai" ]]
  then
    OUTNAME=$(echo "$SAMPLE.$REGIONNAME.cpg.methyl.tsv")
    CPGISLANDFILE="${snakemake_config[CPGISLANDLOC]}/$REGIONNAME.cpg.positions"
    # extract chromosome slash region. if only chromosome don't specify region start end, just use chrom.
    samtools mpileup -r $CHROM --positions $CPGISLANDFILE --ff SUPPLEMENTARY,UNMAP,SECONDARY,QCFAIL,DUP -q 1 -Q 1 -M --no-output-ins --no-output-ins --no-output-del --no-output-del --no-output-ends -o $SAMPLE.$REGIONNAME.pileup.tsv $BAMFILE

    # process the 4th column that contains the data string
    # give the python or rust program a simple input so no expensive parsing is necessary. just a column of strings.

    paste <(cut -f 1,2 --output-delimiter='.' $SAMPLE.$REGIONNAME.pileup.tsv) <(cut -f 5 $SAMPLE.$REGIONNAME.pileup.tsv | ${snakemake_config[RUSTEXE]} T) | sort -k 1,1 -k 2,2n > $OUTNAME
	
	
    rm $SAMPLE.$REGIONNAME.pileup.tsv
	
	
  else
    echo "skipping $BAMFILE because no index file exists."
  fi
}

echo -e "chr.pos\t$SAMPLE" > $OUTPUTNAME

if [ $PARALLEL -eq 1 ]
  then
      for i in $(seq 0 $(( ${#REGIONNAMES[@]} - 1)) )
      do
         echo "${REGIONNAMES[$i]} ${CHROMS[$i]} $SAMPLE"
         extractChromosome ${REGIONNAMES[$i]} ${CHROMS[$i]} $SAMPLE $INPUTBAM &
      done
      wait
      for i in $(seq 0 $(( ${#REGIONNAMES[@]} - 1)) )
      do
         cat $SAMPLE.${REGIONNAMES[$i]}.cpg.methyl.tsv >> $OUTPUTNAME
         rm $SAMPLE.${REGIONNAMES[$i]}.cpg.methyl.tsv
      done
  else
      for i in $( seq 0 $(( ${#REGIONNAMES[@]} - 1 )) )
      do
        echo "${REGIONNAMES[$i]} ${CHROMS[$i]} $SAMPLE"
        extractChromosome ${REGIONNAMES[$i]} ${CHROMS[$i]} $SAMPLE $INPUTBAM
        cat $SAMPLE.${REGIONNAMES[$i]}.cpg.methyl.tsv >> $OUTPUTNAME
        rm $SAMPLE.${REGIONNAMES[$i]}.cpg.methyl.tsv
      done
fi

exit 0