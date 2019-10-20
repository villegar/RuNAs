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
GENOME4STAR_FILENAMES = filenames2(GENOME4STAR.keys(),".gz")
PATTERN_R1 = "{library}_1.fastq.gz".format(library=LIBS)
PATTERN_R2 = "{library}_2.fastq.gz".format(library=LIBS)


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
#		expand("4.STAR"),
#		#expand("4.STAR/{library}_Aligned.sortedByCoord.out.bam", library=LIBS)
#		#expand("4.STAR/{library}_{star_file}",library=LIBS,
#		#	star_file=["Aligned.sortedByCoord.out.bam","Aligned.sortedByCoord.out.bam.bai","Log.final.out","Log.out","Log.progress.out","SJ.out.tab","Unmapped.out.mate1","Unmapped.out.mate2"])
#		#expand("4.STAR/STAR.report.html", library = LIBS)
		"done.txt"
rule reads:	
	input:
		reads = READS + "/{library}_{replicate}.fastq.gz",
		r1    = READS + "/{library}_1.fastq.gz",
		r2    = READS + "/{library}_2.fastq.gz"
#		reads = expand(READS + "/{library}_{replicate}.fastq.gz", library=LIBS, replicate=[1, 2]),
#		r1    = expand(READS + "/{library}_1.fastq.gz", library=LIBS),
#		r2    = expand(READS + "/{library}_2.fastq.gz", library=LIBS)
#	output:
#		directory(READS)
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
	message:
		"Trimming reads"
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
		zip = "3.QC.TRIMMED/{library}_{direction}_{mode}_fastqc.zip"
	shell:
                "fastqc -q -o 3.QC.TRIMMED -t {threads} {input}"

rule download_genome:
	output:
		genome_files = expand("GENOME/{genome_file}", genome_file = GENOME4STAR_FILENAMES)
	#	, step = rules.fastqc_trimmed.output.step
	run:
		for link_index in sorted(GENOME4STAR.keys()):
            		shell("wget -q {link} -O GENOME/{file}".format(link=GENOME4STAR[link_index], file=link_index))
			shell("gunzip GENOME/{file}".format(file=link_index))

rule genome_index:
	input:
		genome_files = rules.download_genome.output.genome_files
	output:
		dir = directory("GENOME_INDEX")
	shell:
		"mkdir -p {output.dir} && STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {output} --genomeFastaFiles {input.genome_files[0]}  --sjdbGTFfile {input.genome_files[1]} --sjdbOverhang 50"
		
rule star:
	input:
		genome = rules.genome_index.output.dir,
		r1 = expand("2.TRIMMED/{library}_forward_paired.fastq.gz", library=LIBS),
		r2 = expand("2.TRIMMED/{library}_reverse_paired.fastq.gz", library=LIBS)
	output:
		directory("4.STAR")
		, step = "10"
	#	"4.STAR/{library}_Aligned.sortedByCoord.out.bam"
	params:
		out_prefix = expand("4.STAR/${library}", library=LIBS)
	shell:
		"STAR --runThreadN {threads} --genomeDir {input.genome} --readFilesIn {input.r1} {input.r2} --readFilesCommand gunzip -c --outFilterIntronMotifs RemoveNoncanonical --outFileNamePrefix {params.out_prefix} --outSAMtype BAM SortedByCoordinate --outReadsUnmapped  Fastx"

rule finish:
	input: 
		rules.star.output.step
#		rules.star.output,
#		rules.fastqc_raw.output,
#		rules.fastqc_trimmed.output
	output:
		"done.txt"
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


