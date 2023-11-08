rule index_test_sample:
    input:
        test_case = "results/{SAMPLE}.cpg.methyl.tsv",
        reference_index = ancient("results/meow.reference.index.tsv")
    output:
        "results/{SAMPLE}.cpg.methyl.indexed.tsv"
    threads: 1
    shell:
        """
        LC_ALL=C
        echo -e "index\t{wildcards.SAMPLE}" > {output}
        join <(tail -n+2 {input.reference_index} | sort -k2,2) <(tail -n+2 {input.test_case} | sort -k1,1) -t $'\t' -e 'NA' -a1 -1 2 -2 1 -o "1.1 2.2" >> {output}
        """

rule merge_test_sample:
    input:
        test_case="results/{SAMPLE}.cpg.methyl.indexed.tsv",
        reference_index = ancient("results/meow.reference.index.tsv"),
        reference_database = ancient("results/meow.reference.annotated.tsv")
    output:
        temp("results/{SAMPLE}.annotated.db.tsv")
    threads: 1
    shell:
        """
        LC_ALL=C
        paste <(cut -f1 {input.reference_index}) <(cat {input.reference_database}) > temp_index_{wildcards.SAMPLE}.tsv
        numFields=$( awk '{{print NF;exit}}' temp_index_{wildcards.SAMPLE}.tsv)
        mergeString1=$(for i in $(seq 2 $(echo "$numFields" | bc)); do echo "1.$i"; done | tr '\n' ' ')
        mergeString2="2.2"
        fullMergeString=$(echo "0 $mergeString1$mergeString2" | xargs)
        header="index $(head -n1 {input.reference_database})\t{wildcards.SAMPLE}"
        echo $header | tr ' ' '\t' > {output}
        join <(tail -n+2 temp_index_{wildcards.SAMPLE}.tsv | sort -k1,1) <(tail -n+2 {input.test_case} | sort -k1,1) -t $'\t' -e 'NA' -a1 -o "$fullMergeString" >> {output}
        rm temp_index_{wildcards.SAMPLE}.tsv
        """