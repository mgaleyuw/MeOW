rule generate_pooled_ttest_rmd:
    input:
        "results/{SAMPLE}.annotated.db.tsv"
    output:
        "results/{SAMPLE}_meow_results/{SAMPLE}.pooled_t.dmr.Rmd"
    params:
        clr="FALSE",
        minPositions=50,
        bootstraps=50
    threads: 1
    shell:
        """
        mkdir -p results/{wildcards.SAMPLE}_meow_results
        cat workflow/scripts/dmr_run.Rmd > {output}
        echo -e "\n" >> {output}
        echo "# {wildcards.SAMPLE} MeOW Results" >> {output}
        echo "## Pooled T-test\n\n" >> {output}
        cat workflow/scripts/dmr_pooled.Rmd >> {output}
        echo -e "\n" >> {output}
        cat workflow/scripts/dmr_embedPlots.Rmd >> {output}
        #cp {output} workflow/scripts/rmd_scripts/{wildcards.SAMPLE}.pooled_t.dmr.Rmd
        """

rule generate_paired_ttest_rmd:
    input:
        "results/{SAMPLE}.annotated.db.tsv"
    output:
        "results/{SAMPLE}_meow_results/{SAMPLE}.paired_t.dmr.Rmd"
    params:
        clr="FALSE",
        minPositions=50,
        bootstraps=50,
        control_metadata=config["control_metadata"],
        mdf=config["metadata"]
    threads: 1
    shell:
        """
        mkdir -p results/{wildcards.SAMPLE}_meow_results
        cat workflow/scripts/dmr_run.Rmd > {output}
        echo -e "\n" >> {output}
        echo "# {wildcards.SAMPLE} MeOW Results" >> {output}
        echo "## Paired T-test\n\n" >> {output}
        sed 's/paired=FALSE/paired=TRUE/g' workflow/scripts/dmr_pooled.Rmd >> {output}
        echo -e "\n" >> {output}
        cat workflow/scripts/dmr_embedPlots.Rmd >> {output}
        mkdir -p workflow/results/{wildcards.SAMPLE}_meow_results
        cp {output} workflow/results/{wildcards.SAMPLE}_meow_results/{wildcards.SAMPLE}.paired_t.dmr.Rmd
        cp {input} workflow/results/
        mkdir -p workflow/config
        cp {params.control_metadata} workflow/{params.control_metadata}
        cp {params.mdf} workflow/{params.mdf}
        """

rule generate_beta_ttest_rmd:
    input:
        "results/{SAMPLE}.annotated.db.tsv"
    output:
        "results/{SAMPLE}_meow_results/{SAMPLE}.beta.dmr.Rmd"
    params:
        minPositions=50,
        bootstraps=50,
        control_metadata=config["control_metadata"],
        mdf=config["metadata"]
    threads: 1
    shell:
        """
        mkdir -p results/{wildcards.SAMPLE}_meow_results
        cat workflow/scripts/dmr_run.Rmd > {output}
        echo -e "\n" >> {output}
        echo "# {wildcards.SAMPLE} MeOW Results" >> {output}
        echo "## Beta regression\n\n" >> {output}
        cat workflow/scripts/dmr_beta.Rmd >> {output}
        echo -e "\n" >> {output}
        cat workflow/scripts/dmr_embedPlots.Rmd >> {output}
        mkdir -p workflow/results/{wildcards.SAMPLE}_meow_results
        cp {output} workflow/results/{wildcards.SAMPLE}_meow_results/{wildcards.SAMPLE}.beta.dmr.Rmd
        cp {input} workflow/results/
        mkdir -p workflow/config
        cp {params.control_metadata} workflow/{params.control_metadata}
        cp {params.mdf} workflow/{params.mdf}
        """

rule run_test:
    input:
        cpg_file="results/{SAMPLE}.annotated.db.tsv",
        rmd_doc="results/{SAMPLE}_meow_results/{SAMPLE}.{TEST}.dmr.Rmd"
    output:
        "results/{SAMPLE}_meow_results/{SAMPLE}.{TEST}.dmr.html"
    conda: "envs/tooldmr.yaml"
    params:
        control_metadata=config["control_metadata"],
        mdf=config["metadata"],
        minPositions=50,
        bootstraps=50,
        file_output_path="results/{SAMPLE}_meow_results"
    script:
        "results/{wildcards.SAMPLE}_meow_results/{wildcards.SAMPLE}.{wildcards.TEST}.dmr.Rmd"
