if "restrict-regions" in config["processing"]:

    rule compose_regions:
        input:
            config["processing"]["restrict-regions"],
        output:
            "results/called/{contig}.regions.bed",
        conda:
            "../envs/bedops.yaml"
        shell:
            "bedextract {wildcards.contig} {input} > {output}"


rule call_variants:
    input:
        bam=get_sample_bams,
        bai = (get_bai
               if not config['processing']['bqsr']
               else []),
        ref=get_genome_fun,
        idx=config["local_genome_copy"]["path_to_genome"] + ".dict" if config["local_genome_copy"]["path_to_genome"] != "" else "resources/genome.dict",
        known=config["local_genome_copy"]["known_variants"],
        tbi=config["local_genome_copy"]["known_variants"] + ".tbi",
        regions=(
            "results/called/{contig}.regions.bed"
            if config["processing"].get("restrict-regions")
            else []
        ),
    output:
        gvcf=protected("results/called/{sample}.{contig}.g.vcf.gz"),
    log:
        "logs/gatk/haplotypecaller/{sample}.{contig}.log",
    params:
        extra=get_call_variants_params,
    wrapper:
        "0.59.0/bio/gatk/haplotypecaller"

rule create_strelka_contigs:
    input:
        genome_dict=rules.genome_dict.output
    output:
        bed="results/strelka_beds/{contig}.bed.gz",
        tbi="results/strelka_beds/{contig}.bed.gz.tbi"
    shell:
        """
        echo "#chr    start   end" > create_strelka_contigs.{wildcards.contig}.tmp
        grep  -P SN:{wildcards.contig}"\t" {input.genome_dict} | cut -f 2,3 | sed 's+LN:+1\t+' | sed 's+SN:++' >> create_strelka_contigs.{wildcards.contig}.tmp
        bgzip create_strelka_contigs.{wildcards.contig}.tmp 
        mv create_strelka_contigs.{wildcards.contig}.tmp.gz {output.bed}
        tabix {output.bed}
        """

rule call_strelka:
    input:
        bam=get_sample_bams,
        bai= (get_bai
                if not config['processing']['bqsr']
                else[]),
        ref=get_genome_fun,
        regions="results/strelka_beds/{contig}.bed.gz"
    output:
        vcf="results/called/{sample}_{contig}/results/variants/variants.vcf.gz"
    params:
        strelka=config['local_soft']['strelka']
    conda:"strelka.yml"
    shell:
        """
        {params.strelka} --bam {input.bam} --referenceFasta \
        {input.ref} --runDir results/called/{wildcards.sample}_{wildcards.contig} --callRegions {input.regions} 
        results/called/{wildcards.sample}_{wildcards.contig}/runWorkflow.py -m local --quiet
        """

rule combine_calls:
    input:
        ref=get_genome_fun,
        gvcfs=expand(
            "results/called/{sample}.{{contig}}.g.vcf.gz", sample=samples.index
        )
    output:
        gvcf="results/called/all.{contig}.g.vcf.gz",
    log:
        "logs/gatk/combinegvcfs.{contig}.log",
    wrapper:
        "0.74.0/bio/gatk/combinegvcfs"


rule genotype_variants:
    input:
        ref=get_genome_fun,
        gvcf="results/called/all.{contig}.g.vcf.gz",
    output:
        vcf=temp("results/genotyped/all.{contig}.vcf.gz"),
    params:
        extra=config["params"]["gatk"]["GenotypeGVCFs"],
    log:
        "logs/gatk/genotypegvcfs.{contig}.log",
    wrapper:
        "0.74.0/bio/gatk/genotypegvcfs"


rule merge_variants:
    input:
        vcfs=lambda w: expand(
            "results/genotyped/all.{contig}.vcf.gz", contig=get_contigs()
        ),
    output:
        vcf="results/genotyped/all.vcf.gz",
    log:
        "logs/picard/merge-genotyped.log",
    wrapper:
        "0.74.0/bio/picard/mergevcfs"


rule merge_single_sample_contigs_after_strelka:
    input:
        vcfs=expand(rules.call_strelka.output,sample="{sample}", contig=get_contigs())
    output:
        "results/called/{sample}.vcf.gz",
    log:
        "logs/picard/{sample}.merge-genotyped.log",
    wrapper:
        "0.74.0/bio/picard/mergevcfs"

rule merge_samples_after_strelka:
    input:
        expand(rules.merge_single_sample_contigs_after_strelka.output, sample=samples.index)
    output:
        vcf="results/genotyped/all.strelka.vcf.gz",
    run:
        from subprocess import call
        if len(samples.index) > 1:
            call("bcftools merge  {input}  |  bcftools norm -m - -O z -o  {output}".format(input=input,
            output=output), shell = True)
            call("tabix {output}".format(output=output),shell=True)
        else:
            call("mv {input} {output} ".format(input=input,
            output=output),shell=True)
            call("tabix {output}".format(output=output),shell=True)
