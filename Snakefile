####### Libraries #######
import glob
import os
#from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider

####### Util functions #######
def filenames(path,prefix,suffix):
	filenames_path = glob.glob(os.path.join(path,prefix) + "*" + suffix)
	names = []
	for file in filenames_path:
		names.append(os.path.basename(file).split(suffix)[0])
	return sorted(names)

def filenames2(paths,suffix):
        names = []
        for file in paths:
                names.append(os.path.basename(file).split(suffix)[0])
        return sorted(names)

def which(file):
        for path in os.environ["PATH"].split(os.pathsep):
                if os.path.exists(os.path.join(path, file)):
                        return path
        return None

####### Global variables #######
#READS = "/gpfs/scratch/Classes/stat736/p53reads"
PREFIX = "SRR"
#EXTENSION = "fastq.gz"
READS = "/gpfs/scratch/Classes/stat736/walnutreads"
EXTENSION = "fastq"
SUFFIX = "_1." + EXTENSION
CPUS_FASTQC = 3
CPUS_PHIX = 15
CPUS_TRIMMING = 5
CPUS_STAR = 20
CPUS_ARIA = 16
CPUS_KRAKEN = 20
CPUS_RNA = 20
LIBS = filenames(READS,PREFIX,SUFFIX)
#LIBS = ["SRR2121770"]
#FTP = FTPRemoteProvider()
ADAPTER = which("trimmomatic")
GENOME4STAR = {
	"GRCm38.primary_assembly.genome.fa.gz" : "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M22/GRCm38.primary_assembly.genome.fa.gz",
	"gencode.vM22.annotation.gtf.gz" : "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M22/gencode.vM22.annotation.gtf.gz"
}
GENOME4STAR_FILENAMES = filenames2(GENOME4STAR.keys(),".gz")

GENOME4PHIX = {
	"PhiX" : "ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/PhiX/Illumina/RTA/PhiX_Illumina_RTA.tar.gz"
}
KRAKEN_DB = {
	"minikraken_20171019_8GB" : "https://ccb.jhu.edu/software/kraken/dl/minikraken_20171019_8GB.tgz"
}
RRNA = {
	"txid9606.fasta" : "https://raw.githubusercontent.com/villegar/RuNAs/v2/txid9606.fasta"
}
rRNA_FILES = list(RRNA.keys())
#print(rRNA_FILES)
#shell("mkdir -p GENOME")
#for file, link in GENOME4STAR.items():
#	shell("curl -sS -L {link} -o GENOME/{file}")
#	shell("gunzip GENOME/{file}")

####### Rules #######
rule all:
	input:
		expand("1.QC.RAW/{library}_{replicate}_fastqc.{format}", library=LIBS, replicate=[1, 2], format=["html","zip"]),
		expand("2.TRIMMED/{library}_{direction}_{mode}.fastq.gz",
                        library=LIBS, direction=["forward","reverse"], mode=["paired","unpaired"]),
		expand("3.QC.TRIMMED/{library}_{direction}_{mode}_fastqc.{format}", 
			library=LIBS, direction=["forward","reverse"], mode=["paired","unpaired"], format=["html","zip"]),
#		#expand("GENOME/{genome_file}", genome_file = GENOME4STAR_FILENAMES),
#		#expand("GENOME_INDEX"),
#		expand("4.STAR", library=LIBS),
#		#expand("4.STAR/{library}_Aligned.sortedByCoord.out.bam", library=LIBS)
		expand("4.STAR/{library}_{star_file}", library=LIBS,
			star_file=["Aligned.sortedByCoord.out.bam","Unmapped.out.mate1","Unmapped.out.mate2"]),
##			star_file=["Aligned.sortedByCoord.out.bam","Aligned.sortedByCoord.out.bam.bai","Log.final.out","Log.out","Log.progress.out","SJ.out.tab","Unmapped.out.mate1","Unmapped.out.mate2"])
		expand("5.PHIX/{library}.sam", library=LIBS),
		expand("6.MICROBIAL/{library}.{format}", library=LIBS, format=["out","tsv"]),
		expand("7.rRNA/{library}.rna.{format}", library=LIBS, format=["bam","sam","out"])

rule reads:	
	input:
		reads = READS + "/{library}_{replicate}." + EXTENSION,
		r1    = READS + "/{library}_1." + EXTENSION,
		r2    = READS + "/{library}_2." + EXTENSION
	message:
		"Gathering reads"

rule fastqc_raw:
	input:
		reads = rules.reads.input.reads
		#, step  = rules.reads.output.step
#		reads = READS + "/{library}_{replicate}.fastq.gz"
	output:	
		html = "1.QC.RAW/{library}_{replicate}_fastqc.html",
		zip  = "1.QC.RAW/{library}_{replicate}_fastqc.zip"
	message:
		"FastQC on raw data"
	threads:
		CPUS_FASTQC
	shell:
		"fastqc -q -o 1.QC.RAW -t {threads} {input}"

rule trim_reads:
	input:
		adapter = os.path.join(ADAPTER,"../share/trimmomatic/adapters"),
#		r1      = READS + "/{library}_1.fastq.gz",
#		r2	= READS + "/{library}_2.fastq.gz"
		r1 = rules.reads.input.r1,
		r2 = rules.reads.input.r2
	output:
		forward_paired   = "2.TRIMMED/{library}_forward_paired.fastq.gz",
		forward_unpaired = "2.TRIMMED/{library}_forward_unpaired.fastq.gz",
		reverse_paired   = "2.TRIMMED/{library}_reverse_paired.fastq.gz",
		reverse_unpaired = "2.TRIMMED/{library}_reverse_unpaired.fastq.gz"
	message:
		"Trimming reads"
	threads:
		CPUS_TRIMMING
	shell:
		"trimmomatic PE -threads {threads} {input.r1} {input.r2} {output.forward_paired} {output.forward_unpaired} {output.reverse_paired} {output.reverse_unpaired} ILLUMINACLIP:{input.adapter}/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36"

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
	threads:
		CPUS_FASTQC
	shell:
                "fastqc -q -o 3.QC.TRIMMED -t {threads} {input}"

#rule download_genome:
##	input:
##		FTP.remote(expand("{link}",link=GENOME4STAR.values()), keep_local = True, immediate_close=True)
#	output:
#		genome_files = expand("GENOME/{genome_file}", genome_file = GENOME4STAR_FILENAMES)
#	run:
##		shell("mv {input} GENOME")
#		#for link_index in sorted(GENOME4STAR.keys()):
##		for file, link in GENOME4STAR.items():
##			shell("curl -sS -L {link} -o GENOME/{file}")
##			shell("gunzip GENOME/{file}")	
##			shell("bash download_files.sh {link} && mv {link} GENOME".format(link=GENOME4STAR[link_index]))
##			shell("mv {input} GENOME")
##            		#shell("wget -q -O - {link} | gunzip -c > GENOME/{file}".format(link=GENOME4STAR[link_index], file = filenames2(link_index,".gz")))
##			shell("wget -q {link} -O GENOME/{file} && gunzip GENOME/{file}".format(link=GENOME4STAR[link_index], file=link_index))
##			shell("gunzip GENOME/{file}".format(file=link_index))

rule genome_index:
	input:
		genome_files = expand("GENOME/{genome_file}", genome_file = GENOME4STAR_FILENAMES)
#		genome_files = rules.download_genome.output.genome_files
	output:
		dir = directory("GENOME_INDEX")
	threads: 
		40
	shell:
		"mkdir -p {output.dir} && STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {output} --genomeFastaFiles {input.genome_files[0]}  --sjdbGTFfile {input.genome_files[1]} --sjdbOverhang 50"
		
rule star:
	input:
		genome = rules.genome_index.output.dir,
#		r1 = "2.TRIMMED/{library}_forward_paired.fastq.gz",
#		r2 = "2.TRIMMED/{library}_reverse_paired.fastq.gz"
		r1 = rules.trim_reads.output.forward_paired,
		r2 = rules.trim_reads.output.reverse_paired
	output:
#		"4.STAR/{library}_{star_file}"
		unmapped_m81 = "4.STAR/{library}_Unmapped.out.mate1",
		unmapped_m82 = "4.STAR/{library}_Unmapped.out.mate2",
		aligned_bam  = "4.STAR/{library}_Aligned.sortedByCoord.out.bam"
#		directory("4.STAR")
	message:
		"STAR alignment"
	params:
		prefix = "4.STAR/{library}_"
	threads:
		CPUS_STAR
	shell:
		"STAR --runThreadN {threads} --genomeDir {input.genome} --readFilesIn {input.r1} {input.r2} --readFilesCommand gunzip -c --outFilterIntronMotifs RemoveNoncanonical --outFileNamePrefix {params.prefix} --outSAMtype BAM SortedByCoordinate --outReadsUnmapped  Fastx"

rule phiX_genome:
	output:
		genome = directory(expand("{genome_phix}", genome_phix = GENOME4PHIX.keys()))
	message:
		"Downloading PhiX genome reference"
	run:
		for link_index in sorted(GENOME4PHIX.keys()):
			shell("wget -q -O - {link} | tar -xz".format(link=GENOME4PHIX[link_index]))

rule phiX_contamination:
	input:
		genome = rules.phiX_genome.output.genome,
#		genome = expand("{genome_phix}/Illumina/RTA/Sequence/Bowtie2Index/genome", genome_phix = GENOME4PHIX.keys()),
#		genome	= expand(rules.phiX_genome.output.genome + "Illumina/RTA/Sequence/Bowtie2Index/genome"),
#		genome 	= rules.phiX_genome.output.genome + ["Illumina/RTA/Sequence/Bowtie2Index/genome"],
		r1 	= "2.TRIMMED/{library}_forward_paired.fastq.gz",
                r2 	= "2.TRIMMED/{library}_reverse_paired.fastq.gz"
	message:
		"PhiX contamination analysis"
	output:
		sam = "5.PHIX/{library}.sam"
	threads:
		CPUS_PHIX
	shell:
		"bowtie2 -p {threads} -x {input.genome}/Illumina/RTA/Sequence/Bowtie2Index/genome -1 {input.r1} -2 {input.r2} -S {output}"	

rule kraken_db:
	output:
		kraken_db = directory(expand({"db"}, db=["KRAKEN_DB"]))
	#	kraken_db = directory(expand("{db}", db = KRAKEN_DB.keys()))
	message:
		"Downloading Kraken DB"
	threads:
		CPUS_ARIA
	run:
		for link_index in sorted(KRAKEN_DB.keys()):
			shell("mkdir -p KRAKEN_DB")
			shell("aria2c -x {threads} -s {threads} -d KRAKEN_DB {link}".format(link=KRAKEN_DB[link_index],threads=CPUS_ARIA))	
		#	shell("mv {link_index} {output.kraken_db}")
			shell("tar -xz KRAKEN_DB/{link_index}")
		#	shell("wget -q -O - {link} | tar -xz".format(link=KRAKEN_DB[link_index]))

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
	threads:
		CPUS_KRAKEN
	shell:
		"krakenuniq --preload --db {input.kraken_db} --threads {threads} --paired --report-file {output.tsv} --fastq-input {input.unmapped_m81} {input.unmapped_m82} > {output.out}"

rule rRNA_index:
	output:
		#index = directory(expand("{bwa}", bwa=["BWA_INDEX"])),
		fasta = directory(expand("{bwa}/{file}", bwa=["BWA_INDEX"],file=rRNA_FILES))
	message:
		"Create rRNA index"
	run:
		for link_index in sorted(RRNA.keys()):
			shell("mkdir -p {output.index}")
			shell("wget -q {link}".format(link=RRNA[link_index]))
			shell("mv {link_index} {output.index}")
#			shell("wget -q {link} && mv {link_index} {output.index}".format(link=RRNA[link_index]))
			shell("bwa index {output.index}/{link_index}")
	
rule rRNA_contamination:
	input:
		r1 = rules.trim_reads.output.forward_paired,
		r2 = rules.trim_reads.output.reverse_paired,
		index = rules.rRNA_index.output.fasta
	output:
		sam = "7.rRNA/{library}.rna.sam",
		bam = "7.rRNA/{library}.rna.bam",
		out = "7.rRNA/{library}.rna.out"
	message:
		"Running BWA and Samtools to find rRNA contamination"
	threads:
		CPUS_RNA
	run:
		shell("bwa mem -t {threads} {input.index} {input.r1} {input.r2} > {output.sam}")
		shell("samtools view -@ {threads} -bS -o {output.bam} {output.sam}")
		shell("samtools flagstat -@ {threads} {output.bam} > {output.out}")
		
#rule finish:
#	input: 
#		rules.star.output
#		rules.star.output,
#		rules.fastqc_raw.output,
#		rules.fastqc_trimmed.output
#	output:
#		"done.txt"
#	shell:
#		"echo 'Done'"

#rule multiqc:
#	input:
#		"4.STAR"
#		#"4.STAR/{library}_Aligned.sortedByCoord.out.bam"
#	output:
#		"4.STAR/STAR.report.html"
#	shell:
#		"multiqc -n {output} --flat {input}"


