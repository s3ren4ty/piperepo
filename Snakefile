########################################
# Samples
########################################

SAMPLES = [
    "SRR5660030",
    "SRR5660033",
    "SRR5660044",
    "SRR5660045"
]

########################################
# Final Target
########################################

rule all:
    input:
        "PipelineReport.txt"

########################################
# Build Bowtie2 Index
########################################

rule build_index:
    input:
        "resources/genome/hcmv.fasta"
    output:
        "resources/genome/hcmv_index.1.bt2"
    shell:
        """
        bowtie2-build resources/genome/hcmv.fasta resources/genome/hcmv_index
        """

########################################
# Count Reads BEFORE Mapping
########################################

rule count_before:
    input:
        "data/raw/{sample}_1.fastq"
    output:
        "results/counts/{sample}_before.txt"
    shell:
        """
        mkdir -p results/counts
        echo $(( $(wc -l < {input}) / 4 )) > {output}
        """

########################################
# Map Reads & Keep Only Mapped
########################################

rule map_reads:
    input:
        r1="data/raw/{sample}_1.fastq",
        r2="data/raw/{sample}_2.fastq",
        index="resources/genome/hcmv_index.1.bt2"
    output:
        "results/mapped/{sample}.bam"
    shell:
        """
        mkdir -p results/mapped
        bowtie2 -x resources/genome/hcmv_index \
            -1 {input.r1} \
            -2 {input.r2} | \
            samtools view -b -F 4 - > {output}
        """

########################################
# Count Reads AFTER Mapping
########################################

rule count_after:
    input:
        "results/mapped/{sample}.bam"
    output:
        "results/counts/{sample}_after.txt"
    shell:
        """
        echo $(( $(samtools view -c {input}) / 2 )) > {output}
        """

########################################
# Convert BAM to Paired FASTQ
########################################

rule bam_to_fastq:
    input:
        "results/mapped/{sample}.bam"
    output:
        r1="results/filtered_fastq/{sample}_1.fastq",
        r2="results/filtered_fastq/{sample}_2.fastq"
    shell:
        """
        mkdir -p results/filtered_fastq
        samtools sort -n {input} -o results/mapped/{wildcards.sample}_sorted.bam
        samtools fastq \
            -1 {output.r1} \
            -2 {output.r2} \
            -0 /dev/null -s /dev/null -n \
            results/mapped/{wildcards.sample}_sorted.bam
        """

########################################
# Assemble with SPAdes (k=99)
########################################

rule assemble_spades:
    input:
        r1="results/filtered_fastq/{sample}_1.fastq",
        r2="results/filtered_fastq/{sample}_2.fastq"
    output:
        "results/assembly/{sample}/contigs.fasta"
    shell:
        """
        mkdir -p results/assembly/{wildcards.sample}
        spades.py \
            -1 {input.r1} \
            -2 {input.r2} \
            -k 99 \
            -o results/assembly/{wildcards.sample}
        """

########################################
# Assembly Statistics (>1000 bp)
########################################

rule assembly_stats:
    input:
        "results/assembly/{sample}/contigs.fasta"
    output:
        "results/assembly_stats/{sample}_assembly.txt"
    run:
        from Bio import SeqIO
        import os

        os.makedirs("results/assembly_stats", exist_ok=True)

        count = 0
        total_length = 0

        for record in SeqIO.parse(input[0], "fasta"):
            length = len(record.seq)
            if length > 1000:
                count += 1
                total_length += length

        with open(output[0], "w") as out:
            out.write(f"{count}\t{total_length}")

########################################
# Build Final Report
########################################

rule build_report:
    input:
        before=expand("results/counts/{sample}_before.txt", sample=SAMPLES),
        after=expand("results/counts/{sample}_after.txt", sample=SAMPLES),
        assembly=expand("results/assembly_stats/{sample}_assembly.txt", sample=SAMPLES)
    output:
        "PipelineReport.txt"
    run:
        with open(output[0], "w") as out:
            for s in SAMPLES:
                with open(f"results/counts/{s}_before.txt") as b:
                    before = b.read().strip()

                with open(f"results/counts/{s}_after.txt") as a:
                    after = a.read().strip()

                with open(f"results/assembly_stats/{s}_assembly.txt") as f:
                    count, total = f.read().strip().split("\t")

                out.write(
                    f"Sample {s} had {before} read pairs before and {after} read pairs after Bowtie2 filtering.\n"
                )

                out.write(
                    f"In the assembly of sample {s}, there are {count} contigs > 1000 bp and {total} total bp.\n\n"
                )
