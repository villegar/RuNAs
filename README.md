
# R<sub>u</sub>NA-s<sub>eq</sub> (&#5809;<sub>&#5794;</sub>&#5822;&#5800;&#5835;) */roonas/* Pipeline

![Rule Graph](images/rule-graph.png?raw=true "Rule Graph")
## Requirements
-	BEDtools 2.29+
-	Bowtie2 2.3.5+
-	BWA (alignment via Burrows-Wheeler transformation) 0.7.17+
-	FastQC 0.11.8+
-	FeatureCounts 1.6.4+
-	KrakenUniq 0.5.8+
-	MultiQC 1.7+
-	Python 3.6+ (using [Ana](https://anaconda.org)([mini](https://docs.conda.io/en/latest/miniconda.html))conda)
-	Samtools 1.9+
-	Snakemake 5.7+
-	SRA Toolkit 2.9.6+
-	STAR 2.7.3+
-	Trimmomatic 0.39+

## Setup
```bash
git clone https://github.com/villegar/runas
cd runas
conda env create -f environment.yml
conda activate RuNAs or source activate RuNAs
python download.genome.py human-genome.json
```

## Execution
### Single node
```bash
snakemake -j CPUS \ # maximum number of CPUs available to Snakemake
	  --configfile config.json # configuration file
```

### Multi-node
```bash
snakemake -j JOBS  \ # maximum number of simultaneous jobs to spawn
	  --configfile config.json # configuration file
          --latency-wait 1000 \ # files latency in seconds
          --cluster-config cluster.json \ # cluster configuration file
          --cluster "sbatch --job-name={cluster.name} 
                            --nodes={cluster.nodes} 
                            --ntasks-per-node={cluster.ntasks} 
                            --output={cluster.log} 
                            --partition={cluster.partition} 
                            --time={cluster.time}"
```
#### Alternatively
```bash
bash run_cluster &> log &
```

#### Cluster configuration (cluster.json)
```bash
{
    "__default__" :
    {
        "time" : "1-00:00:00",
        "nodes" : 1,
        "partition" : "compute",
	"ntasks": "{threads}",
	"name": "RuNAs-{rule}",
	"log": "RuNAS-{rule}-%J.log"
    }
}
```

#### Pipeline configuration (config.json)
```bash
{
    "genome4phiX":
    {
        "PhiX": 
            "ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/PhiX/Illumina/RTA/PhiX_Illumina_RTA.tar.gz"
    },

    "genome4star":
    {
        "GRCm38.primary_assembly.genome.fa.gz": 
            "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M22/GRCm38.primary_assembly.genome.fa.gz",
        "gencode.vM22.annotation.gtf.gz": 
            "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M22/gencode.vM22.annotation.gtf.gz"
    },

    "krakenDB":
    {
        "minikraken_20171019_8GB.tgz": 
            "https://ccb.jhu.edu/software/kraken/dl/minikraken_20171019_8GB.tgz"
    },
    "_comment" :
        "Below READS section shows the configuration for a directory containing reads in the format:
	/{PATH}/SRR{LIBRARY}_1.fastq",
    "reads":
    {
        "extension": "fastq.gz",
        "path": "/gpfs/scratch/Classes/stat736/p53reads",
        "prefix": "SRR"
    },

    "rRNAref":
    {
        "txid9606.fasta": 
            "https://raw.githubusercontent.com/villegar/RuNAs/v2/txid9606.fasta"
    }
}
```
