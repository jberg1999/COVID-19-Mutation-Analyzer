# Covid-19-Mutation-Analyzer
A website for exploratory data analysis of Covid-19 genomes

Link to website (https://jberg1999.shinyapps.io/COVID19MutationAnalyzer/?_ga=2.7382962.1609937937.1618995377-1954151250.1614141124)

Hey yall,

Thanks for taking the time to look over this project!. The website, and in turn the repository, is split up into two working parts. 
The first part is the pipeline for aligning and processing sequences. Unfortunately, due to the size limits on GitHub, I am unable 
to post the actual sequence file required to run the processing locally. However it can be recreated by downloading sequences from 
NCBI Virus 
(https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Severe%20acute%20respiratory%20syndrome%20coronavirus%202,%20taxid:2697049)
in fasta format and merging them with metadata available from the same site. Below is a breakdown of each file and what it does.


Backend

processMaster.py: takes sequence file, splits it, and creates threads of process.py to handle each part. Needs a sequence file

Process.py: takes a split sequence file, aligns the sequences, processes the alignments for mutations, and stores the sequence 
and mutation data into files. Needs refseq.fasta and a file created by processMaster.


Front end

App.R: main file that runs the website. Needs Functions.R, app_mut, app_seq, and refseq.fasta

Functions.R: contains functions required for App.R to run

app_mut.csv: a csv file containing mutations and their metadata that is read into app.R

app_seq.csv: a csv file containing sequences and their metadata that is read into app.R


Other files

refseq.fasta: a fasta file containing the NCBI reference genome for SARS-CoV-2, NC_045512, which was collected from Wuhan, China 
in December 2019. It is the oldest known covid sequence. Both process and app.R require this to run

walkthrough poster: a PowerPoint file containing an overview of the project as well as a quick demo of how to use the website.
