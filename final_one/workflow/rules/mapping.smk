rule trim_reads_se:
    input:
        unpack(get_fastq),
    output:
        temp("results/trimmed/{sample}-{unit}.fastq.gz"),
    params:
        **config["params"]["trimmomatic"]["se"],
        extra="",
    log:
        "logs/trimmomatic/{sample}-{unit}.log",
    wrapper:
        "0.74.0/bio/trimmomatic/se"


rule trim_reads_pe:
    input:
        unpack(get_fastq),
    output:
        r1=temp("results/trimmed/{sample}-{unit}.1.fastq.gz"),
        r2=temp("results/trimmed/{sample}-{unit}.2.fastq.gz"),
        r1_unpaired=temp("results/trimmed/{sample}-{unit}.1.unpaired.fastq.gz"),
        r2_unpaired=temp("results/trimmed/{sample}-{unit}.2.unpaired.fastq.gz"),
        trimlog="results/trimmed/{sample}-{unit}.trimlog.txt",
    params:
        **config["params"]["trimmomatic"]["pe"],
        extra=lambda w, output: "-trimlog {}".format(output.trimlog),
    log:
        "logs/trimmomatic/{sample}-{unit}.log",
    wrapper:
        "0.74.0/bio/trimmomatic/pe"


rule fastp_se:
    input:
        unpack(get_fastq)
    output:
        trimmed="results/fast[/{sample}.fastq.gz",
        html="report/se/{sample}.html",
    log:
        "logs/fastp/se/{sample}.log"
    shell:
        "fastp --in1 {input.r1} --out1 {output.out1} -h {output.html} 2> {log}"



rule fastp_pe:
    input:
        unpack(get_fastq)
    output:
        out1 = temp("results/fastp/{sample}-{unit}.1.fastq.gz"),
        out2 = temp("results/fastp/{sample}-{unit}.2.fastq.gz"),
        # Unpaired reads separately
        unpaired1=temp("results/fastp/{sample}-{unit}.1.unpaired.fastq.gz"),
        unpaired2=temp("results/fastp/{sample}-{unit}.2.unpaired.fastq.gz"),
        # or in a single file
#        unpaired="trimmed/pe/{sample}.singletons.fastq",
        html="report/pe/{sample}-{unit}.html",
    log:
        "logs/fastp/pe/{sample}-{unit}.log"
    shell:
        "fastp --in1 {input.r1} --in2 {input.r2} --out1 {output.out1} --out2 {output.out2} --unpaired1 {output.unpaired1} --unpaired2 {output.unpaired2} -h {output.html} 2> {log}"

rule map_reads_with_minimap:
    input:
        reads=get_trimmed_reads_fastp if config['processing']['fastp'] else get_trimmed_reads,
        idx=rules.minimap_index.output
    output:
        temp("results/mapped/{sample}-{unit}.sorted.minimap.bam")
    log:
        "logs/minimap/{sample}-{unit}.log"
    shell:
        """
         minimap2 -Y  -ax sr  -R '@RG\\tID:~{wildcards.sample}\\tSM:~{wildcards.sample}\\tPL:BGI\\tLB:lib1\\tPU:unit1' {input.idx} {input.reads}  |  samtools sort -O sam -o results/mapped/test.sam
         samtools view -b -h results/mapped/test.sam -o {output} 
        """

rule map_reads:
    input:
        reads=get_trimmed_reads,
        idx=rules.bwa_index.output,
    output:
        temp("results/mapped/{sample}-{unit}.sorted.bam"),
    log:
        "logs/bwa_mem/{sample}-{unit}.log",
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
        extra=get_read_group,
        sort="samtools",
        sort_order="coordinate",
    threads: 8
    wrapper:
        "0.74.0/bio/bwa/mem"


rule mark_duplicates:
    input:
        "results/mapped/{sample}-{unit}.sorted.minimap.bam" if config["processing"]["minimap"] else "results/mapped/{sample}-{unit}.sorted.bam"
    output:
        bam=temp("results/dedup/{sample}-{unit}.bam"),
        metrics="results/qc/dedup/{sample}-{unit}.metrics.txt",
    log:
        "logs/picard/dedup/{sample}-{unit}.log",
    params:
        config["params"]["picard"]["MarkDuplicates"],
    wrapper:
        "0.74.0/bio/picard/markduplicates"


rule recalibrate_base_qualities:
    input:
        bam=get_recal_input(),
        bai=get_recal_input(bai=True),
        ref=get_genome_fun,
        dict=config["local_genome_copy"]["path_to_genome"] + ".dict" if config["local_genome_copy"]["path_to_genome"] != "" else "resources/genome.dict",
        known=config["local_genome_copy"]["known_variants"],
        known_idx=config["local_genome_copy"]["known_variants"] + ".tbi",
    output:
        recal_table="results/recal/{sample}-{unit}.grp",
    log:
        "logs/gatk/bqsr/{sample}-{unit}.log",
    params:
        extra=get_regions_param() + config["params"]["gatk"]["BaseRecalibrator"],
    resources:
        mem_mb=1024,
    wrapper:
        "0.74.0/bio/gatk/baserecalibrator"


rule apply_base_quality_recalibration:
    input:
        bam=get_recal_input(),
        bai=get_recal_input(bai=True),
        ref=get_genome_fun,
        dict=config["local_genome_copy"]["path_to_genome"] + ".dict" if config["local_genome_copy"]["path_to_genome"] != "" else "resources/genome.dict",
        recal_table="results/recal/{sample}-{unit}.grp",
    output:
        bam=protected("results/recal/{sample}-{unit}.bam"),
    log:
        "logs/gatk/apply-bqsr/{sample}-{unit}.log",
    params:
        extra=get_regions_param(),
    resources:
        mem_mb=1024,
    wrapper:
        "0.74.0/bio/gatk/applybqsr"


rule samtools_index_dedup:
    input:
        "results/dedup/{sample}-{unit}.bam",
    output:
        "results/dedup/{sample}-{unit}.bam.bai",
    log:
        "logs/samtools/index_mdup/{sample}-{unit}.log",
    wrapper:
        "0.74.0/bio/samtools/index"


