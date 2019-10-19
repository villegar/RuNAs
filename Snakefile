#LIBS=["SRR2121770","SRR2121771","SRR2121774"]

####### Libraries #######
import glob
import os

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

READS = "/gpfs/scratch/Classes/stat736/p53reads"
PREFIX = "SRR"
SUFFIX = "_1.fastq.gz"
LIBS = filenames(READS,PREFIX,SUFFIX)
ADAPTER = which("trimmomatic")
GENOME4STAR = {
	"GRCm38.primary_assembly.genome.fa.gz" : "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M22/GRCm38.primary_assembly.genome.fa.gz",
	"gencode.vM22.annotation.gtf.gz" : "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M22/gencode.vM22.annotation.gtf.gz"
}

####### Rules #######

rule all:
	input:
		expand("1.QC.RAW/{library}_{replicate}_fastq.html", library=LIBS, replicate=[1, 2]),
		expand("1.QC.RAW/{library}_{replicate}_fastq.zip", library=LIBS, replicate=[1, 2]),
		expand("2.TRIMMED/trimm_{library}_forward_paired.fastq.gz", library=LIBS),
		expand("2.TRIMMED/trimm_{library}_forward_unpaired.fastq.gz", library=LIBS),
		expand("2.TRIMMED/trimm_{library}_reverse_paired.fastq.gz", library=LIBS),
		expand("2.TRIMMED/trimm_{library}_reverse_unpaired.fastq.gz", library=LIBS),
		expand("3.QC.TRIMMED/{library}_{direction}_{mode}_fastq.html", 
			library=LIBS, direction=["forward","reverse"], mode=["paired","unpaired"]),
                expand("3.QC.TRIMMED/{library}_{direction}_{mode}_fastq.zip", 
			library=LIBS, direction=["forward","reverse"], mode=["paired","unpaired"]),
		expand("4.STAR/{library}_Aligned.sortedByCoord.out.bam",library=LIBS)
		#expand("4.STAR/{library}_{star_file}",library=LIBS,
		#	star_file=["Aligned.sortedByCoord.out.bam","Aligned.sortedByCoord.out.bam.bai","Log.final.out","Log.out","Log.progress.out","SJ.out.tab","Unmapped.out.mate1","Unmapped.out.mate2"])

rule fastqc_raw:
	input:
		r1 = "reads/{library}_1.fastq.gz",
		r2 = "reads/{library}_2.fastq.gz"
	output:	
		"1.QC.RAW/{library}_{replicate}_fastq.html",
		"1.QC.RAW/{library}_{replicate}_fastq.zip"
	shell:
		"fastqc -o 1.QC.RAW -t {threads} {input}"

rule trimm_reads:
	input:
		adapter = os.path.join(ADAPTER,"../share/trimmomatic/adapters"),
		r1 = "reads/{library}_1.fastq.gz",
                r2 = "reads/{library}_2.fastq.gz"
	output:
		forward_paired = "2.TRIMMED/trimm_{library}_forward_paired.fastq.gz",
		forward_unpaired = "2.TRIMMED/trimm_{library}_forward_unpaired.fastq.gz",
		reverse_paired = "2.TRIMMED/trimm_{library}_reverse_paired.fastq.gz",
                reverse_unpaired = "2.TRIMMED/trimm_{library}_reverse_unpaired.fastq.gz"
	shell:
		"trimmomatic PE -threads {threads} {input.r1} {input.r2} {output.forward_paired} {output.forward_unpaired} {output.reverse_paired} {output.reverse_unpaired} ILLUMINACLIP:{input.adapter}/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36"

rule fastqc_trimmed:
	input:
		"2.TRIMMED/trimm_{library}_{direction}_{mode}.fastq.gz"
		#forward_paired = "2.TRIMMED/trimm_{library}_forward_paired.fastq.gz",
                #forward_unpaired = "2.TRIMMED/trimm_{library}_forward_unpaired.fastq.gz",
                #reverse_paired = "2.TRIMMED/trimm_{library}_reverse_paired.fastq.gz",
                #reverse_unpaired = "2.TRIMMED/trimm_{library}_reverse_unpaired.fastq.gz"
	output:
                "3.QC.TRIMMED/{library}_{direction}_{mode}_fastq.html",
                "3.QC.TRIMMED/{library}_{direction}_{mode}_fastq.zip"
	shell:
                "fastqc -o 3.QC.TRIMMED -t {threads} {input}"
rule download_genome:
	output:
		expand("GENOME/{file}", file = GENOME4STAR.keys())
	run:
		for link_index in sorted(GENOME4STAR.keys()):
            		shell("wget {link} -O GENOME/{file}".format(link=GENOME4STAR[link_index], file=link_index))
			shell("gunzip GENOME/{file}".format(file=link_index))

GENOME4STAR_FILENAMES = filenames2(GENOME4STAR.keys(),".gz")
rule genome_index:
	input:
		genome_files = expand("GENOME/{file}", file = GENOME4STAR_FILENAMES)
	output:
		directory("GENOME_INDEX")
	shell:
		"STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {output} --genomeFastaFiles {input.genome_files[0]}  --sjdbGTFfile {input.genome_files[1]} --sjdbOverhang 50"
		
rule star:
	input:
		genome = "GENOME",
		r1 = "2.TRIMMED/trimm_{library}_forward_paired.fastq.gz",
		r2 = "2.TRIMMED/trimm_{library}_reverse_paired.fastq.gz"
	output:
		"4.STAR/{library}_Aligned.sortedByCoord.out.bam"
	#	prefix = "{library}_Aligned.sortedByCoord.out.bam"
	shell:
		"STAR --runThreadN {threads} --genomeDir {input.genome} --readFilesIn {input.r1} {input.r2} --readFilesCommand gunzip -c --outFilterIntronMotifs RemoveNoncanonical --outFileNamePrefix 4.STAR --outSAMtype BAM SortedByCoordinate --outReadsUnmapped  Fastx"


