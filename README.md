
#in order to successfully run this pipeline, we will need to have all of the following installed.

#Snakemake 
#Bowtie2
#SAMtools
#SPAdes
#BLAST+
#Python
#Biopython 
#Seqio 

#we will also need raw paired end reads in fastq file format. I was able to acquire these samples using the ‘prefetch’ command and then entering in my sample number (EX:prefetch SRR5660030) and then used fasterq-dump --split-files to split the files. I also first made a directory that contains all of my raw data for easy access.

#It is also important to note that the current sample data in my repo are just a small subset of the input reads (10,000), as i shortened them (referencing this source https://genomicislands.wordpress.com/2013/11/20/subsampling-large-fastq-files/#:~:text=Subsampling%20large%20fastq%20files%20If%20you're%20doing,an%20idea%20of%20how%20the%20analysis%20works.)


#if you want to run your own paired end data, please put it in the raw data file.

#in order to run the code, input “snakemake --cores 4” into your terminal

#sources:

#-# https://en.wikipedia.org/wiki/SAMtools
#- https://snakemake.readthedocs.io/en/stable/tutorial/basics.html#step-1-mapping-reads
#- https://submit.mit.edu/submit-users-guide/tutorials/tutorial_7.html
#- https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.blastx_application_options/#:~:text=The%20blastx%20application%20translates%20a,sequences%20or%20a%20protein%20database.
#- https://ifb-elixirfr.github.io/IFB-FAIR-bioinfo-training/assets/pdf/Session2020/session_03/03_workflow.pdf
#- https://www.youtube.com/watch?v=2SbeficAisg
#- 
# pipeproj
