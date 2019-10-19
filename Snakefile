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


####### Rules #######

rule all:
	input:
		expand("1.QC.RAW/{library}_{replicate}_fastq.html", library=LIBS, replicate=[1, 2]),
		expand("1.QC.RAW/{library}_{replicate}_fastq.zip", library=LIBS, replicate=[1, 2]),
		expand("2.TRIMMED/trimm_{library}_forward_paired.fastq.gz", library=LIBS),
		expand("2.TRIMMED/trimm_{library}_forward_unpaired.fastq.gz", library=LIBS),
		expand("2.TRIMMED/trimm_{library}_reverse_paired.fastq.gz", library=LIBS),
		expand("2.TRIMMED/trimm_{library}_reverse_unpaired.fastq.gz", library=LIBS)

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
                r1 = "2.TRIMMED/{library}_1.fastq.gz",
                r2 = "2.TRIMMED/{library}_2.fastq.gz"
        output:
                "1.QC.RAW/{library}_{replicate}_fastq.html",
                "1.QC.RAW/{library}_{replicate}_fastq.zip"
        shell:
                "fastqc -o 1.QC.RAW -t {threads} {input}"
