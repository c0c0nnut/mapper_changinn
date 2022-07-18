rule get_genome:
    output:
        "resources/genome.fasta",
    log:
        "logs/get-genome.log",
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    cache: True
    wrapper:
        "0.74.0/bio/reference/ensembl-sequence"


checkpoint genome_faidx:
    input:
        get_genome_fun,
    output:
        config["local_genome_copy"]["path_to_genome_fai"] if config["local_genome_copy"]["path_to_genome_fai"] != "" else "resources/genome.fasta.fai",
    log:
        "logs/genome-faidx.log",
    cache: True
    wrapper:
        "0.74.0/bio/samtools/faidx"


rule genome_dict:
    input:
        get_genome_fun,
    output:
        config["local_genome_copy"]["path_to_genome"] + ".dict" if config["local_genome_copy"]["path_to_genome"] != "" else "resources/genome.dict"
    log:
        "logs/samtools/create_dict.log",
    conda:
        "../envs/samtools.yaml"
    cache: True
    shell:
        "samtools dict {input} > {output} 2> {log} "


rule get_known_variation:
    input:
        # use fai to annotate contig lengths for GATK BQSR
        fai=get_fai,
    output:
        vcf="resources/variation.vcf.gz",
    log:
        "logs/get-known-variants.log",
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        type="all",
    cache: True
    wrapper:
        "0.74.0/bio/reference/ensembl-variation"


# rule remove_iupac_codes:
#     input:
#         "resources/variation.vcf.gz",
#     output:
#         "resources/variation.noiupac.vcf.gz",
#     log:
#         "logs/fix-iupac-alleles.log",
#     conda:
#         "../envs/rbt.yaml"
#     cache: True
#     shell:
#         "rbt vcf-fix-iupac-alleles < {input} | bcftools view -Oz > {output}"
#
#
# rule tabix_known_variants:
#     input:
#         "resources/variation.noiupac.vcf.gz",
#     output:
#         "resources/variation.noiupac.vcf.gz.tbi",
#     log:
#         "logs/tabix/variation.log",
#     params:
#         "-p vcf",
#     cache: True
#     wrapper:
#         "0.74.0/bio/tabix"


rule bwa_index:
    input:
        get_genome_fun,
    output:
        multiext(config["local_genome_copy"]["path_to_genome"] if config["local_genome_copy"]["path_to_bwa_index"] != "" else "resources/genome.fasta" , ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa_index.log",
    resources:
        mem_mb=369000,
    cache: True
    wrapper:
        "0.74.0/bio/bwa/index"


rule minimap_index:
    input:
        get_genome_fun,
    output:
        multiext(config["local_genome_copy"]["path_to_genome"] if config["local_genome_copy"]["path_to_minimap_index"] != "" else "resources/genome.fasta" , ".mni"),
    shell:
        "minimap2 -d {output} {input}  "

rule get_vep_cache:
    output:
        directory("resources/vep/cache"),
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    log:
        "logs/vep/cache.log",
    wrapper:
        "0.74.0/bio/vep/cache"


rule get_vep_plugins:
    output:
        directory("resources/vep/plugins"),
    log:
        "logs/vep/plugins.log",
    params:
        release=config["ref"]["release"],
    wrapper:
        "0.74.0/bio/vep/plugins"
