
# RuNAs (&#5809;&#5794;&#5822;&#5800;&#5835;): RN<sub>u</sub>As<sub>eq</sub> */roonas/* Pipeline

![Rule Graph](https://raw.github.com/villegar/RuNAs/v2/images/rule-graph.svg?raw=true "Rule Graph")
## Requirements
-	Python 3.6+ (using [Ana](https://anaconda.org)([mini](https://docs.conda.io/en/latest/miniconda.html))conda)

## Setup
```bash
git clone https://github.com/villegar/runas --branch v2
cd runas
conda env create -f environment.yml
python download.genome.py human-genome.json
snakemake -j CPUS 
```
