# inDrops

## Installation

Installer avec conda l'environnement `single_cell_snakemake.yaml` ainsi que l'environnement `single_cell.yaml`

```bash
conda env create -f single_cell_snakemake.yaml
```

## YAML file configuration

Aide pour régler les problèmes liés au fichier YAML pour plus de détails voir [la page github correspondante](https://github.com/indrops/indrops).

Indiquer le nom du projet et l'emplacement des fichiers de sortie.

```yaml
project_name : "SingleCell_ICM"
project_dir : "/home/adrien.dufour/PROTECT/debug_output/"
```

```yaml
cores :
  default : 3
  quantify_barcodes : 8
```

Nombre de coeurs pour les différentes étapes, augmenter le nombre de coeurs de `quantify_barcodes` en cas d'erreur samtools

```bash
[E::hts_open_format] Failed to open file "scRNA_project/sample1/quant_dir/bcIKBG.genomic.sorted.bam" : Too many open files
samtools merge: fail to open "scRNA_project/sample1/quant_dir/bcIKBG.genomic.sorted.bam": Too many open files
```

Adapter le nom des fichiers fasta à la nomenclature de la version désirée (v1, v2, v3).

```yaml
sequencing_runs :
  - name : "Run1"
    version : 'v1'
    dir : "/home/adrien.dufour/PROTECT/debug_data/"
    fastq_path : "{library_prefix}_{split_affix}_{read}.fastq.gz"
    split_affixes : ["S3", "S4"]
    libraries :
      - {library_name: "AR005", library_prefix: "AR005"}
      - {library_name: "AR004", library_prefix: "AR004"}
```

Indiquer l'emplacement des fichiers d'index en `.annotated`

Indiquer l'emplacement du `/bin/` de l'environnement `single_cell_snakemake.yaml` pour Bowtie, Samtools, Rsem et Java

Indiquer l'emplacement du `/bin/` de l'environnement `single_cell.yaml` pour Python

```yaml
paths :
  bowtie_index : "/home/adrien.dufour/PROTECT/DOWNLOAD_DIR/Mus_musculus.GRCm38.85.annotated"
  bowtie_dir : "/opt/miniconda3/bin/"
  python_dir : "/opt/miniconda3/envs/indrops/bin/"
  samtools_dir : "/opt/miniconda3/bin/"
  rsem_dir : "/opt/miniconda3/bin/"
  java_dir : "/opt/miniconda3/bin/"
```

Paramètres additionnels

```yaml
parameters : # OPTIONAL PARAMETERS
  umi_quantification_arguments:
    m : 10 #Ignore reads with more than M alignments, after filtering on distance from transcript end.
    u : 1 #Ignore counts from UMI that should be split among more than U genes.
    d : 400 #Maximal distance from transcript end, NOT INCLUDING THE POLYA TAIL
    split-ambigs: False #If umi is assigned to m genes, add 1/m to each gene's count (instead of 1)
    min_non_polyA: 0 #Require reads to align to this much non-polyA sequence. (Set to 0 to disable filtering on this parameter.)
  output_arguments:
    output_alignment_to_bam: True
    output_unaligned_reads_to_other_fastq: False
    output_oversequencing_metrics: True
    output_umifm_calculation_metrics: True
    output_quant_metrics: True
  bowtie_arguments:
    m : 200
    n : 1
    l : 15
    e : 50
  trimmomatic_arguments:
    LEADING: "28"
    SLIDINGWINDOW: "4:20"
    MINLEN: "16"
```

## Snakemake file configuration

Indiquer l'emplacement du fichier d'environnement `single_cell.yaml`

Indiquer l'emplacement du fichier yaml de configuration `project.yaml`

```python
shell.executable("/bin/bash")
import itertools

conda: "/home/adrien.dufour/NeuroDev_ADD/Envs/single_cell.yaml"
configfile: "/home/adrien.dufour/NeuroDev_ADD/SingleCell/indrops-master/project.yaml"
#shell.prefix("conda activate indrops")

YAMLPATH = "/home/adrien.dufour/NeuroDev_ADD/SingleCell/indrops-master/project.yaml"
```

Modifier en fonction du naming de vos fichiers

```python
RUN_LIBRARY = []
FASTQ = []
LIBRARY = []
WORKERS = range(config['cores']['default'])
WORKER = ['1', '2', '3']
READS = ['R1', 'R2']
SPLIT = []
RUN = []

for each in config['sequencing_runs']:
    RUN = each['name']
    dir_lib = each['dir']
    SPLIT = each['split_affixes']
    RUN_LIBRARY.append((RUN, each['library_name']))
    LIBRARY.append(each['library_name'])
```

## Usage

Ce placer dans le dossier contenant le fichier `indrops.py`

```bash
snakemake --cores 12 --use-conda --conda-prefix "/opt/miniconda3" --latency-wait 25
```

Commande pour le cluster de Strasbourg

```bash
snakemake  \
--configfile /b/home/inci/mokhtari/SingleCell/strasbourg/indrops-master/project-Copy7.yaml \
--cluster-config /b/home/inci/mokhtari/SingleCell/Script/cluster_config_sc.yaml \
--snakefile /b/home/inci/mokhtari/SingleCell/strasbourg/indrops-master/Snakefile_7 \
--jobs 1 \
--cluster "sbatch --verbose \
-A ${ACCOUNT} \
-p ${PARTITION} \
--time={cluster.walltime} \
--mem={cluster.mem_gb}G \
--cpus-per-task={cluster.cpus} \
--output={cluster.stdout} \
--error={cluster.stderr}" \
--use-conda \
--conda-prefix "/b/home/inci/mokhtari/.conda/envs/singlecell/bin/" \
--printshellcmds \
--latency-wait 60 \
--keep-going
```

## Troubleshooting

Si votre environnement de travail ne supporte pas les fichiers fifo modifiez la ligne 1231 du fichier `indrops.py`

```python
filtered_dir = os.path.dirname("/mnt/ram")
```
