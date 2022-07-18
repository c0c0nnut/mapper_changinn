
rule create_pileup_call:
    input:
        index=rules.index_bam.output,
        bam=rules.sort_sam.output,
        genome="resources/genome.fasta"
    output:
        "results/calls/{sample}.vcf"
    log:
        "logs/samtools/{sample}.mpileup.log"
    params:
        call_params=config["params"]["call_params"]
    shell:
        """
        samtools mpileup -g -f {input.genome} {input.bam} | bcftools call {params.call_params} - > {output}
        """
