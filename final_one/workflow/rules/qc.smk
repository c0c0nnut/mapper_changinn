rule fastqc:
    input:
        unpack(get_fastq),
    output:
        html="results/qc/fastqc/{sample}-{unit}.html",
        zip="results/qc/fastqc/{sample}-{unit}.zip",
    log:
        "logs/fastqc/{sample}-{unit}.log",
    wrapper:
        "0.74.0/bio/fastqc"


rule samtools_stats:
    input:
        get_sample_bams,
    output:
        "results/qc/samtools-stats/{sample}-{unit}.txt",
    log:
        "logs/samtools-stats/{sample}-{unit}.log",
    wrapper:
        "0.74.0/bio/samtools/stats"


rule multiqc:
    input:
        expand(
            [
                "results/qc/samtools-stats/{u.sample}-{u.unit}.txt",
                "results/qc/fastqc/{u.sample}-{u.unit}.zip",
                "results/qc/dedup/{u.sample}-{u.unit}.metrics.txt",
            ],
            u=units.itertuples(),
        ),
    output:
        report(
            "results/qc/multiqc.html",
            caption="../report/multiqc.rst",
            category="Quality control",
        ),
    log:
        "logs/multiqc.log",
    run:
        from os import path
        from subprocess import call

        output_dir = path.dirname(output[0])
        output_name = path.basename(output[0])

        call("multiqc --force -o {output_dir} -n {output_name} {input_data} > {log} 2>&1".format(output_dir = output_dir,
        output_name = output_name, input_data = input[0], log = log ), shell = True )
    #wrapper:
    #    "0.74.0/bio/multiqc"
