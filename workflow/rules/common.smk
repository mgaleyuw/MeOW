import pandas as pd
import glob

samples=pd.read_table(config["samples"], sep="\t", index_col=0, header=0)
samples["output_name"]=samples.index+"."+samples["Type"]+".allChr.tsv"

chrs = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]

def get_input_bam(wildcards):
    return samples.loc[wildcards.SAMPLE, "Bamfile"]


def get_all_ttest_targets(wildcards):
    #samples["output_name"]=samples.index+"."+samples["Type"]+".allChr.tsv"
    f = open(config["targetfile"], "r")
    samples = f.read().split("\n")
    f.close()
    one_column_file="results/{sample}_meow_results/{sample}.paired_t.dmr.Rmd"
    return [one_column_file.format(sample=x) for x in samples]

def get_all_beta_targets(wildcards):
    f = open(config["targetfile"], "r")
    samples = f.read().split("\n")
    f.close()
    one_column_file="results/{sample}_meow_results/{sample}.beta.dmr.Rmd"
    return [one_column_file.format(sample=x) for x in samples]

def get_one_beta_target(wildcards):
    return "results/{sample}_meow_results/{sample}.beta.dmr.html".format(sample=config["single_target"])

def get_one_ttest_target(wildcards):
    return "results/{sample}_meow_results/{sample}.beta.dmr.html".format(sample=config["single_target"])

def get_all_controls(wildcards):
    bamfiles = samples[samples["Type"]=="Control"].index.tolist()
    return ["results/{samplename}.cpg.methyl.tsv".format(samplename=x) for x in bamfiles]
    

def get_target_type(wildcards):
    return samples.loc[wildcards.SAMPLE, "type"]

def get_rust_exe(methylation):
    match methylation:
        case "5mc":
            return "workflow/scripts/./methyl_pileup"
        case "5hmc5mc_ignore":
            return "workflow/scripts/./methyl_pileup_5hmc_5mcBlind"
        case "5hmc5mc_bayes":
            return "workflow/scripts/./methyl_pileup_5hmc_5mc"
        case _:
            print("option {} not recognized, exiting.".format(methylation))
            exit

rule build_column:
    input:
        chr_methyls = expand("results/{SAMPLE}.{region}.cpg.methyl.tsv", region=chrs, allow_missing=True)
    output: 
        OUTPUTNAME=temp("results/{SAMPLE}.cpg.methyl.tsv")
    threads: 1
    conda: "../envs/toolparse.yaml"
    params:
        CpGFile=config["DEFAULTREGION"]
    shell:
        """
        LC_ALL=C

        echo -e "chr.pos\t{wildcards.SAMPLE}" > {output.OUTPUTNAME}
        BASHLIST=( {input.chr_methyls} )    
        for FILE in ${{BASHLIST[@]}}
        do
            cat $FILE >> {output.OUTPUTNAME}
        done
        """

rule extract_chromosome:
    input:
        bam_file = get_input_bam,
        CPGISLANDFILE="".join([config["CPGISLANDLOC"],"/{REGIONNAME}.cpg.positions"])
    output:
        chr_methyl = temp("results/{SAMPLE}.{REGIONNAME}.cpg.methyl.tsv")
    threads: 1
    conda: "../envs/toolparse.yaml"
    params:
        excl_flags="SUPPLEMENTARY,UNMAP,SECONDARY,QCFAIL,DUP",
        rust_exe=get_rust_exe(config["methylation_type"])
    shell:
        """
        LC_ALL=C
        samtools mpileup -r {wildcards.REGIONNAME} --positions {input.CPGISLANDFILE} --ff {params.excl_flags} -q 1 -Q 1 -M --no-output-ins --no-output-ins --no-output-del --no-output-del --no-output-ends -o scratch/{wildcards.SAMPLE}.{wildcards.REGIONNAME}.pileup.tsv {input.bam_file}
        pwd
        echo "{params.rust_exe}"
        paste <(cut -f 1,2 --output-delimiter='.' scratch/{wildcards.SAMPLE}.{wildcards.REGIONNAME}.pileup.tsv) <(cut -f 5 scratch/{wildcards.SAMPLE}.{wildcards.REGIONNAME}.pileup.tsv | {params.rust_exe} T) | sort -k 1,1 -k 2,2n > {output.chr_methyl}
        rm scratch/{wildcards.SAMPLE}.{wildcards.REGIONNAME}.pileup.tsv
        """
    
