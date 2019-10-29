####### Libraries #######
import glob
import os

####### Util functions #######
def extractFilenames(fullnames,suffix):
        names = []
        for file in fullnames:
            names.append(os.path.basename(file).split(suffix)[0])
        return sorted(names)

def findLibraries(path,prefix,suffix):
	filenames_path = glob.glob(os.path.join(path,prefix) + "*" + suffix)
	names = []
	for file in filenames_path:
	    library = os.path.basename(file).split(suffix)[0]
	    if(library not in names):
		    names.append(library)
	return sorted(names)

def which(file):
        for path in os.environ["PATH"].split(os.pathsep):
                if os.path.exists(os.path.join(path, file)):
                        return path
        return None

####### Global variables #######
EXTENSION = config["reads"]["extension"]
PREFIX = config["reads"]["prefix"]
READS = config["reads"]["path"]
FORWARD_READ_ID = config["reads"]["forward_read_id"]
SUFFIX = "_" + FORWARD_READ_ID + "." + EXTENSION
LIBS = findLibraries(READS,PREFIX,SUFFIX)

###### Multithread configuration #####
CPUS_FASTQC = 3
CPUS_PHIX = 15
CPUS_TRIMMING = 5
CPUS_STAR = 20
CPUS_ARIA = 16
CPUS_KRAKEN = 20
CPUS_READCOUNTS = 5
CPUS_RNA = 20

ADAPTER = which("trimmomatic")

####### Reference datasets #######
GENOME4STAR = config["genome4star"]
GENOME4STAR_FILENAMES = extractFilenames(GENOME4STAR.keys(),".gz")
GENOME4PHIX = config["genome4phiX"]
KRAKEN_DB = config["krakenDB"]
KRAKEN_DB_FILENAMES = extractFilenames(KRAKEN_DB.keys(),".tgz")
RRNA = config["rRNAref"]
rRNA_FILES = list(RRNA.keys())

####### Rules #######
rule all:
	input:
		expand("1.QC.RAW/{library}_{end}_fastqc.{format}", library=LIBS, end=[1, 2], format=["html","zip"]),
		#expand("2.TRIMMED/{library}_{direction}_{mode}.fastq.gz",
                #        library=LIBS, direction=["forward","reverse"], mode=["paired","unpaired"]),
		expand("3.QC.TRIMMED/{library}_{direction}_{mode}_fastqc.{format}", 
			library=LIBS, direction=["forward","reverse"], mode=["paired","unpaired"], format=["html","zip"]),
		#expand("4.STAR/{library}_{star_file}", library=LIBS,
		#	star_file=["Aligned.sortedByCoord.out.bam","Unmapped.out.mate1","Unmapped.out.mate2"]),
		"readCounts.txt",
		expand("5.PHIX/{library}.sam", library=LIBS),
		expand("6.MICROBIAL/{library}.{format}", library=LIBS, format=["out","tsv"]),
		#expand("7.rRNA/{library}.rRNA.{format}", library=LIBS, format=["bam","sam","out"]),
		#expand("8.rRNA.FREE.READS/{library}_{end}.fastq", library=LIBS, end=[1, 2]),
		expand("9.QC.rRNA.FREE.READS/{library}_{end}_fastqc.{format}", library=LIBS, end=[1, 2], format=["html","zip"])
	
	output:
		logs 	= directory("0.LOGS"),
		reports	= directory("10.MULTIQC")
	run:
		shell("multiqc -o {output.reports} -n 1.Report_FastQC_Raw.html -d 1.QC.RAW")
		shell("multiqc -o {output.reports} -n 2.Report_Trimming.html -d 2.TRIMMED")
		shell("multiqc -o {output.reports} -n 3.Report_FastQC_Trimmed.html -d 3.QC.TRIMMED")
		shell("multiqc -o {output.reports} -n 4.Report_STAR.html -d 4.STAR")
		shell("multiqc -o {output.reports} -n 5.Report_PhiX.html -d 5.PHIX")
		shell("multiqc -o {output.reports} -n 6.Report_Microbial.html -d 6.MICROBIAL")
		shell("multiqc -o {output.reports} -n 7.Report_rRNA.html -d 7.rRNA")
		shell("multiqc -o {output.reports} -n 8.Report_rRNA_free.html -d 8.rRNA.FREE.READS")
		shell("multiqc -o {output.reports} -n 9.Report_FastQC_rRNA_free.html -d 9.QC.rRNA.FREE.READS")
		shell("mkdir -p {output.logs} && mv *.log {output.logs}")

rule reads:	
	input:
		reads = READS + "/{library}_{end}." + EXTENSION,
		r1    = READS + "/{library}_1." + EXTENSION,
		r2    = READS + "/{library}_2." + EXTENSION
	message:
		"Gathering reads"

rule fastqc_raw:
	input:
		reads = rules.reads.input.reads
	output:	
		html = "1.QC.RAW/{library}_{end}_fastqc.html",
		zip  = "1.QC.RAW/{library}_{end}_fastqc.zip"
	message:
		"FastQC on raw data"
	log:
		"1.QC.RAW/{library}_{end}.log"
	threads:
		CPUS_FASTQC
	shell:
		"fastqc -o 1.QC.RAW -t {threads} {input} 2> {log}"

rule trim_reads:
	input:
		adapter = os.path.join(ADAPTER,"../share/trimmomatic/adapters"),
		r1 = rules.reads.input.r1,
		r2 = rules.reads.input.r2
	output:
		forward_paired   = "2.TRIMMED/{library}_forward_paired.fastq.gz",
		forward_unpaired = "2.TRIMMED/{library}_forward_unpaired.fastq.gz",
		reverse_paired   = "2.TRIMMED/{library}_reverse_paired.fastq.gz",
		reverse_unpaired = "2.TRIMMED/{library}_reverse_unpaired.fastq.gz"
	message:
		"Trimming reads"
	log:
		"2.TRIMMED/{library}.log"
	threads:
		CPUS_TRIMMING
	shell:
		"trimmomatic PE -threads {threads} {input.r1} {input.r2} {output.forward_paired} {output.forward_unpaired} {output.reverse_paired} {output.reverse_unpaired} ILLUMINACLIP:{input.adapter}/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads SLIDINGWINDOW:4:20 TRAILING:3 MINLEN:36 2> {log}"

rule fastqc_trimmed:
	input:
		rules.trim_reads.output.forward_paired,
		rules.trim_reads.output.forward_unpaired,
		rules.trim_reads.output.reverse_paired,
		rules.trim_reads.output.reverse_unpaired
	output:
		html = "3.QC.TRIMMED/{library}_{direction}_{mode}_fastqc.html",
		zip  = "3.QC.TRIMMED/{library}_{direction}_{mode}_fastqc.zip"
	message:
		"FastQC on trimmed data"
	log:
		"3.QC.TRIMMED/{library}_{direction}_{mode}.log"
	threads:
		CPUS_FASTQC
	shell:
                "fastqc -o 3.QC.TRIMMED -t {threads} {input} 2> {log}"

rule genome_index:
	input:
		genome_files = expand("GENOME/{genome_file}", genome_file = GENOME4STAR_FILENAMES)
	output:
		dir = directory("GENOME_INDEX")
	message:
		"Generate genome index for STAR"
	log:
		"GENOME/genome_index.log"
	threads: 
		CPUS_STAR
	shell:
		"mkdir -p {output.dir} && STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {output} --genomeFastaFiles {input.genome_files[0]}  --sjdbGTFfile {input.genome_files[1]} --sjdbOverhang 50 2> {log}"
		
rule star:
	input:
		genome = rules.genome_index.output.dir,
		r1 = rules.trim_reads.output.forward_paired,
		r2 = rules.trim_reads.output.reverse_paired
	output:
		unmapped_m81 = "4.STAR/{library}_Unmapped.out.mate1",
		unmapped_m82 = "4.STAR/{library}_Unmapped.out.mate2",
		aligned_bam  = "4.STAR/{library}_Aligned.sortedByCoord.out.bam"
	message:
		"STAR alignment"
	log:
		"4.STAR/{library}.log"
	params:
		prefix = "4.STAR/{library}_"
	threads:
		CPUS_STAR
	shell:
		"STAR --runThreadN {threads} --genomeDir {input.genome} --readFilesIn {input.r1} {input.r2} --readFilesCommand gunzip -c --outFilterIntronMotifs RemoveNoncanonical --outFileNamePrefix {params.prefix} --outSAMtype BAM SortedByCoordinate --outReadsUnmapped  Fastx 2> {log}"

rule read_counts:
	input:
		aligned = expand(rules.star.output.aligned_bam, library=LIBS),
		genome = rules.genome_index.input.genome_files[1]
	output:
		readCounts = "readCounts.txt"
	log:
		"read_counts.log"
	threads:
		CPUS_READCOUNTS
	shell:
		"featureCounts -a {input.genome} -o {output} -T {threads} {input.aligned} 2> {log}"
		
rule phiX_genome:
	output:
		genome = directory(expand("{genome_phix}", genome_phix = GENOME4PHIX.keys()))
	message:
		"Downloading PhiX genome reference"
	log:
		"phiX_genome.log"
	run:
		for link_index in sorted(GENOME4PHIX.keys()):
			shell("wget -q -O - {link} | tar -xz".format(link=GENOME4PHIX[link_index]))

rule phiX_contamination:
	input:
		genome 	= rules.phiX_genome.output.genome,
		r1	= rules.trim_reads.output.forward_paired,
                r2 	= rules.trim_reads.output.reverse_paired
	message:
		"PhiX contamination analysis"
	output:
		sam = "5.PHIX/{library}.sam"
	log:
		"5.PHIX/{library}.log"
	threads:
		CPUS_PHIX
	shell:
		"bowtie2 -p {threads} -x {input.genome}/Illumina/RTA/Sequence/Bowtie2Index/genome -1 {input.r1} -2 {input.r2} -S {output} 2> {log}"	

rule kraken_db:
	output:
		kraken_db = directory(expand("KRAKEN_DB/{db}", db = KRAKEN_DB_FILENAMES))
	message:
		"Downloading Kraken DB"
	log:
		"KRAKEN_DB/kraken_db.log"
	threads:
		CPUS_ARIA
	run:
		for link_index in sorted(KRAKEN_DB.keys()):
			shell("aria2c -x {threads} -s {threads} -d KRAKEN_DB {link}".format(link=KRAKEN_DB[link_index],threads=CPUS_ARIA))	
			shell("tar -xzf KRAKEN_DB/{link_index} -C KRAKEN_DB")
			shell("build_taxdb {output.kraken_db}/taxonomy/names.dmp {output.kraken_db}/taxonomy/nodes.dmp > {output.kraken_db}/taxDB 2> {log}")

rule microbial_contamination:
	input:
		unmapped_m81 	= rules.star.output.unmapped_m81,
		unmapped_m82 	= rules.star.output.unmapped_m82,
		kraken_db	= rules.kraken_db.output.kraken_db
	output:
		out = "6.MICROBIAL/{library}.out",
		tsv = "6.MICROBIAL/{library}.tsv"
	message:
		"Running KrakenSeq to find microbial contamination"
#	log:
#		"6.MICROBIAL/{library}.log"
	threads:
		CPUS_KRAKEN
	shell:
		"krakenuniq --preload --db {input.kraken_db} --threads {threads} --paired --report-file {output.tsv} --fastq-input {input.unmapped_m81} {input.unmapped_m82} > {output.out}"

rule rRNA_index:
	output:
		fasta = expand("{bwa}/{file}", bwa=["BWA_INDEX"],file=rRNA_FILES)
	message:
		"Create rRNA index"
	log:
		"rRNA_index.log"
	run:
		for link_index in sorted(RRNA.keys()):
			shell("wget -q {link}".format(link=RRNA[link_index]))
			shell("mv {link_index} BWA_INDEX")
			shell("bwa index {output.fasta} 2> {log}")
	
rule rRNA_contamination:
	input:
		r1 = rules.trim_reads.output.forward_paired,
		r2 = rules.trim_reads.output.reverse_paired,
		index = rules.rRNA_index.output.fasta
	output:
		sam = "7.rRNA/{library}.rRNA.sam",
		bam = "7.rRNA/{library}.rRNA.bam",
		unmapped_bam = "7.rRNA/{library}_unmapped.rRNA.bam",
		out = "7.rRNA/{library}.rRNA.out"
	message:
		"Running BWA and Samtools to find rRNA contamination"
	threads:
		CPUS_RNA
	run:
		shell("bwa mem -t {threads} {input.index} {input.r1} {input.r2} > {output.sam}")
		shell("samtools view -@ {threads} -bS -o {output.bam} {output.sam}")
		shell("samtools flagstat -@ {threads} {output.bam} > {output.out}")
		shell("samtools view -@ {threads} -u -f 12 -F 256 {output.bam} > {output.unmapped_bam}")

rule trim_rRNA_contamination:
	input:
		unmapped_bam = rules.rRNA_contamination.output.unmapped_bam
	output:
		r1 = "8.rRNA.FREE.READS/{library}_1.fastq",
		r2 = "8.rRNA.FREE.READS/{library}_2.fastq"
	log:
		"8.rRNA.FREE.READS/{library}.log"
	shell:
		"bedtools bamtofastq -i {input} -fq {output.r1} -fq2 {output.r2} 2> {log}"

rule fastqc_trimmed_rRNA:
	input:
		r1 = rules.trim_rRNA_contamination.output.r1,
		r2 = rules.trim_rRNA_contamination.output.r2
	output:
		html = "9.QC.rRNA.FREE.READS/{library}_{end}_fastqc.html",
		zip  = "9.QC.rRNA.FREE.READS/{library}_{end}_fastqc.zip"
	message:
		"FastQC on reads without rRNA contamination"
	log:
		"9.QC.rRNA.FREE.READS/{library}_{end}.log"
	threads:
		CPUS_FASTQC
	shell:
		"fastqc -o 9.QC.rRNA.FREE.READS -t {threads} {input} 2> {log}"

