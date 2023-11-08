rule merge_controls:
    input:
        control_files = get_all_controls
    output:
        merge_tsv = "results/meow.reference.merged.tsv"
    threads: 1
    params:
        minPositions=config["minimum_control_CpGs"]
    shell:
        """
        ALLFILES=( {input.control_files} )
        PASSFILES=( )
        for file in ${{ALLFILES[@]}}
        do
            POSITIONS=$(awk 'BEGIN{{nas=0}}NR>1{{if($2=="NA"){{nas=nas+1}}}}END{{print NR-nas}}' $file)
            if [[ $POSITIONS -ge {params.minPositions} ]]
            then
                PASSFILES+=( $file )
                echo "$file passes with $POSITIONS cpg positions"
            else
                echo "$file fails with $POSITIONS cpg positions"
            fi
        done
        
        firstFile=${{PASSFILES[0]}}
        for ifile in ${{PASSFILES[@]:1}}
        do
            secondFile=$ifile
            numFields=$( awk '{{print NF;exit}}' $firstFile)
            mergeString1=$(for i in $(seq 2 $(echo "$numFields" | bc)); do echo "1.$i"; done | tr '\n' ' ')
            numFields=$( awk '{{print NF;exit}}' $secondFile)
            mergeString2=$(for i in $(seq 2 $(echo "$numFields" | bc)); do echo "2.$i"; done | tr '\n' ' ')
            fullMergeString=$(echo "0 $mergeString1$mergeString2" | xargs)
            LC_ALL=C
            echo $fullMergeString
            join <(head -n 1 $firstFile) <(head -n 1 $secondFile) -t $'\t' > {output.merge_tsv}
            LC_ALL=C
            join <(tail -n+2 $firstFile | sort -k 1,1) <(tail -n+2 $secondFile | sort -k 1,1) -t $'\t' -e 'NA' -a1 -a2 -o "$fullMergeString" >> {output.merge_tsv}
            cp {output.merge_tsv} temp_FirstFile.tsv
            firstFile=temp_FirstFile.tsv
        done
        rm temp_FirstFile.tsv
        """

rule annotate_controls:
    input:
        ancient("results/meow.reference.merged.tsv")
    output:
        "results/meow.reference.annotated.tsv"
    threads: 1
    conda: "../envs/toolparse.yaml"
    params:
        rustannotate=config["RUSTANNOTATE"],
        namedcpgs=config["NAMEDCPGS"],
        namedgenes=config["NAMEDGENES"]
    shell:
        """
        echo -e "chr.pos\tcpg" > intermediate.names.tsv
        cat {input} | sed 's/\./\t/' |tail -n+2 | sort -k1,1 -k2,2n | cut -f 1,2 | {params.rustannotate} {params.namedcpgs} h | sort -k1,1 >> intermediate.names.tsv
        echo -e "$(head -n 1 {input})\tcpg" > intermediate.cpgs.named.all.tsv
        join <(tail -n+2 {input} | sort -k1,1) <(tail -n+2 intermediate.names.tsv | sort -k1,1) -t $'\t' >> intermediate.cpgs.named.all.tsv
        # produces a two column file 
        #get gene names
        cat {input} | sed 's/\./\t/' |tail -n+2 | sort -k1,1 -k2,2n | cut -f 1,2 | {params.rustannotate} {params.namedgenes} h | sort -k1,1 > intermediate.gene.names.tsv
        echo -e "$(head -n 1 intermediate.cpgs.named.all.tsv)\tgene" > cpgs.named.all.tsv
        join <(tail -n+2 intermediate.cpgs.named.all.tsv | sort -k1,1) <(tail -n+2 intermediate.gene.names.tsv | sort -k1,1) -t $'\t' >> cpgs.named.all.tsv

        #finally fix the indexing
        paste <(cat cpgs.named.all.tsv | cut -f 1 | tr '.' '\t') <(cut -f2- cpgs.named.all.tsv) > {output}

        rm intermediate*tsv
        rm cpgs.named.all.tsv
        """

rule index_controls:
    input:
        ancient("results/meow.reference.annotated.tsv")
    output:
        "results/meow.reference.index.tsv"
    threads: 1
    shell:
        """
        numlines=$(echo "$(cat {input} | wc -l) -1" | bc )
        paste <(seq 0 $numlines) <(cut -f 1,2 {input} | tr '\t' '.') > {output}
        """