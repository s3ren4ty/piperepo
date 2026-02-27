#List of samples that defines the four diff datasets the pipeline will run on
SAMPLES = ["SRR5660030",
    "SRR5660033",
    "SRR5660044",
    "SRR5660045"]

#Rule all: the final goal of the pipeline is to produce our pipeline report.
#snakemake will work backwards from this to figure out what needs to run
rule all:
    input:
        "PipelineReport.txt"

#Building our bowtienindex from our HCMV gemone fasta file as ourn input (Q2)
rule build_index:
    input:
        "resources/genome/hcmv.fasta" 
        #hcmv fasta as input
    output:
        "resources/genome/hcmv_index.1.bt2" 
        #uses bowtie to build multiple index files that share the same prefix. Allows bowtie to map sequencing reads. (literally like an index in the book. allows bowtie to run efficiently)
    shell:
        """
        bowtie2-build resources/genome/hcmv.fasta resources/genome/hcmv_index """ 
        #shell command 

#rule that counts the reads BEFORE we map them out (Q2)
rule count_before:
    input:
        "data/raw/{sample}_1.fastq" 
        #takes in each samples paired end fastq file 
    output:
        "results/counts/{sample}_before.txt" 
        #output file
    shell:
        """
        mkdir -p results/counts
        echo $(( $(wc -l < {input}) / 4 )) > {output} """ 
        #takes each sample and divides it by 4 since there are 4 lines per read in a fastq file and it is only counting the lines. n ot the full reads

#rule that will map out our reads (Q2)
#takes in our paired end reads for each sample as input. as well as its index
rule map_reads:
    input:
        r1="data/raw/{sample}_1.fastq",
        r2="data/raw/{sample}_2.fastq",
        index="resources/genome/hcmv_index.1.bt2"
    output:
        "results/mapped/{sample}.bam" 
        #where our results will go. in a bam file for each sample
        #shell makes a directory for all the mapped reads
        #bowtie maps the paired end reads to the HCMV genome
        #samtools view -b is used to conver the sam files to bam
        # -F 4 removes all unmapped reads. so the BAM outputs only contains reads successfully mapepd to the HCMV index
    shell:
        """
        mkdir -p results/mapped
        bowtie2 -x resources/genome/hcmv_index \
            -1 {input.r1} \
            -2 {input.r2} | \
            samtools view -b -F 4 - > {output}"""

#Rule that will count the total alignments in each bam file (Q2)
#before filtering count and after filtering count
rule count_after:
    input:
        "results/mapped/{sample}.bam" 
        #input is the bam file for each sample
    output:
        "results/counts/{sample}_after.txt" 
        #output is the sample after it counts the total alignments in bam then divides by two because each pair produces 2 alignments
    shell:
        """
        echo $(( $(samtools view -c {input}) / 2 )) > {output}""" 
        #Uses samtools to view the input file and then divide it by two so they pair together and produce 2 alignments after filtering

#converting BAM files to paired end FASTQ reads (Q3)
rule bam_to_fastq:
    input:
        "results/mapped/{sample}.bam" 
        #input file. using sample bam files
    output:
        r1="results/filtered_fastq/{sample}_1.fastq",
        r2="results/filtered_fastq/{sample}_2.fastq" 
        #outputs are our two filtered paired end reads in fastq format
    shell:
      #making new directory for our filtered fastq reads called filtered_fastq
      #using samtools sort -n to sort BAM by read name, which is needed in a paired FASTQ input
        """
        mkdir -p results/filtered_fastq
        samtools sort -n {input} -o results/mapped/{wildcards.sample}_sorted.bam
        samtools fastq \
            -1 {output.r1} \
            -2 {output.r2} \
            -0 /dev/null -s /dev/null -n \
            results/mapped/{wildcards.sample}_sorted.bam""" 
            #samtools fastq is used to convert the filtered BAM back into a paired fastq file. as SPAdes needs a fastq file for an input

#rule that assembles with SPADES!! using k=99 as per the question asks (Q3)
rule assemble_spades:
    input:
        r1="results/filtered_fastq/{sample}_1.fastq", 
        r2="results/filtered_fastq/{sample}_2.fastq" 
        #input files are the paired end filtered fastq files we produced
    output:
        "results/assembly/{sample}/contigs.fasta" 
        #output file
    shell:
        """
        mkdir -p results/assembly/{wildcards.sample}
        spades.py \
            -1 {input.r1} \
            -2 {input.r2} \
            -k 99 \
            -o results/assembly/{wildcards.sample}""" 
            #spades.py -1 and -2 running spades using paired end files. -k 99 sets imer size to 99. -o specifies output directory


#now using python code
#gets our assembly stats (Q4)
rule assembly_stats:
    input:
        "results/assembly/{sample}/contigs.fasta" 
        #input file. using contig file for each sample
    output:
        "results/assembly_stats/{sample}_assembly.txt" 
        #outputting an assembly file for each sample

        #tells snakemake to run python command
    run: 
        from Bio import SeqIO 
        #importing seqio so we can read fasta files and parse sequences one at a time
        import os 
        #importing os module so we can create directories

        os.makedirs("results/assembly_stats", exist_ok=True) 
        #creating our assembly_stats directory. exist_ok=true preventing an error if the folder exists already (so python doesnt crash...)
        #initializing two counters
        count = 0 
        #counter for the number of contigs > 1000 bp
        total_length = 0 
        #sum of the contig lengths

        #input[0] takes the first input file (contigs,fasta file)
        #'fasta' tells SEQIO the file format. before using it t parse the file one sequence at a time
        #each record represnents one contig with attributes like record.id and record. seq
        for record in SeqIO.parse(input[0], "fasta"): 
            length = len(record.seq) 
            #getting length of contig sequence
            if length > 1000: 
              #only taking contigs longer than 1000 bp
                count += 1 
                increase count by 1
                total_length += length 
                #adding that contigs length to the total length

        #writing output to output file defined earlier
        #wriiting two numbers separated by a tab. the number of contigs, and the total length
        with open(output[0], "w") as out:
            out.write(f"{count}\t{total_length}")

#rule that will extract our longest contig (Q5)
#using python now
rule extract_longest_contig:
    input:
        "results/assembly/{sample}/contigs.fasta" 
        #input file contigs.fasta
    output:
        "results/blast/{sample}_longest.fasta" 
        #output file the longest contig found 
    run:
        from Bio import SeqIO #using SEQIO 
        import os #using os again

        os.makedirs("results/blast", exist_ok=True) 
        #making a new directory for BLAST results

        longest = None 
        #creating variable to store the longest contig in it
        max_len = 0 
        #making an variable to track the length of the longest contig seen so far. starting at 0 as it will take in squences with a longer length

        for record in SeqIO.parse(input[0], "fasta"): 
          #reading in the input fasta file parsing it one contig at a time
        #each record is a seqrecord object
            if len(record.seq) > max_len: 
              #for each contig, measure its length and compare it to the current longest contig. if its longer make it the current max length holder
                longest = record 
                #storign the entire contig as the longest
                max_len = len(record.seq) 
                #updating the maximum contig length value so the program can keep track

        SeqIO.write(longest, output[0], "fasta") 
        #writing out longest contig to an output file. taking whatever was in our 'longest; variable

# rule that makes the Betaherpesvirinae BLAST database (Q5)
rule make_betaherpes_db:
    input:
        "resources/blast/betaherpesvirinae.fasta" 
        #taking in our betaherpesvirinae fasta file
    shell:
        """
        makeblastdb -in {input} -dbtype nucl -out resources/blast/betaherpesvirinae_db""" 
        #shell executed code. creating a local BLAST database limited to only betaherpesvirinae as the question asks

#rule that will BLAST our longest contig
rule blast_longest_contig:
    input:
        query="results/blast/{sample}_longest.fasta" 
        #taking in our longest contig fasta file for each sample
    output:
        "results/blast/{sample}_blast.txt" 
        #our output file
    #our shell commands
    #using blastn as our query and database are both nucleotide
    # -max_hsps 1 keeps only best alignment per query subject pair as question asks for while keeping the top 5 hits
    # -outfmt is our custom tab-delimited output. contains all the stuff question 5 requires
    shell:
        """
        blastn \
          -query {input.query} \
          -db resources/blast/betaherpesvirinae_db \
          -max_hsps 1 \
          -num_alignments 5 \
          -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle" \
          > {output}""" 

#rule that builds our final PipelineReport
#using python
#combining all the results from our previous rules into one file!
rule build_report:
    #expand() takes list of samples. using (SAMPLES) to generate all file paths. using our paired end reads like S1 and S2 for each sample
    #taking our before, after, assembly, and blast inputs
    input:
        before=expand("results/counts/{sample}_before.txt", sample=SAMPLES),
        after=expand("results/counts/{sample}_after.txt", sample=SAMPLES),
        assembly=expand("results/assembly_stats/{sample}_assembly.txt", sample=SAMPLES),
        blast=expand("results/blast/{sample}_blast.txt", sample=SAMPLES)

    #our output file
    output:
        "PipelineReport.txt"

    #running python code when rule is execuded
    #looping overeach sample in SAMPLES (our sample list). 
    run:
        with open(output[0], "w") as out:
            for s in SAMPLES:
                
                #taking in read pairs after filtering and storing value in before
                with open(f"results/counts/{s}_before.txt") as b:
                    before = b.read().strip()

                #taking in read pairs after filtering and storing value in after
                with open(f"results/counts/{s}_after.txt") as a:
                    after = a.read().strip()

                #opening assembly stats for sample. reading in content and then splitting line at the tab into two variables. counting the number of contigs > 1000 bp and the total length of those contigs
                with open(f"results/assembly_stats/{s}_assembly.txt") as f:
                    count, total = f.read().strip().split("\t")

                #human readable lines of output
                out.write(
                    f"Sample {s} had {before} read pairs before and {after} read pairs after Bowtie2 filtering.\n")

              
                out.write(
                    f"In the assembly of sample {s}, there are {count} contigs > 1000 bp and {total} total bp.\n\n")
                out.write(f"{s}:\n")
                #writing our blast table headers 
                out.write("sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n")

                #opening the BLAST output for each sample and looping over each line of BLAST output. writing it into the report under the header
                with open(f"results/blast/{s}_blast.txt") as blast:
                    for line in blast:
                        out.write(line)

                #adding double blank lines btween samples just for clarity
                out.write("\n\n")
