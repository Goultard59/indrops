shell.executable("/bin/bash")
import itertools

conda: "/home/adrien.dufour/NeuroDev_ADD/Envs/single_cell.yaml"
configfile: "/home/adrien.dufour/NeuroDev_ADD/SingleCell/indrops-master/project.yaml"
#shell.prefix("conda activate indrops")

YAMLPATH = "/home/adrien.dufour/NeuroDev_ADD/SingleCell/indrops-master/project.yaml"
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
    LIBRARY.append(each['library_name'])

def aggregate_input(wildcards):
    library_quant = [os.path.join(config['project_dir'], wildcards.library, "quant_dir",
                     "worker{i}_".format(i=i) + str(config['cores']['quantify_barcodes']) + ".counts.tsv") \
                    for i in WORKERS]
    return library_quant

rule all:
    input:
        [os.path.join(config['project_dir'], 'fastqc', x.replace('.fastq', '_fastqc.html')) for x in FASTQ if 'R1' in x],
        expand(os.path.join(config['project_dir'], "{library}", "filtered_parts", "{library}_{run}_{split}.fastq"), split=SPLIT, library=LIBRARY, run=RUN),
        expand(os.path.join(config['project_dir'], "{library}", "filtered_parts", "{library}_{run}_{split}.fastq.counts.pickle"), split=SPLIT, library=LIBRARY, run=RUN),
        expand(os.path.join(config['project_dir'], "{library}", "filtered_parts", "{library}_{run}_{split}metrics.yaml"), split=SPLIT, library=LIBRARY, run=RUN),
        expand(os.path.join(config['project_dir'], "{library}", "abundant_barcodes.pickle"), split=SPLIT, library=LIBRARY, run=RUN),
        expand(os.path.join(config['project_dir'], "{library}", "{library}.barcode_abundance_by_barcode.png"), split=SPLIT, library=LIBRARY, run=RUN),
        expand(os.path.join(config['project_dir'], "{library}", "{library}.barcode_abundance.png"), library=LIBRARY),
        expand(os.path.join(config['project_dir'], "{library}", "{library}.filtering_stats.csv"), library=LIBRARY),
        expand(os.path.join(config['project_dir'], "{library}", "filtered_parts", "{library}_{run}_{split}.fastq.sorted.fastq.gz"), split=SPLIT, run=RUN, library=LIBRARY),
        expand(os.path.join(config['project_dir'], "{library}", "filtered_parts", "{library}_{run}_{split}.fastq.sorted.fastq.gz.index.pickle"), split=SPLIT, run=RUN, library=LIBRARY),
        expand(os.path.join(config['project_dir'], "{library}", "quant_dir", "worker{i}_" + str(config['cores']['quantify_barcodes']) + ".counts.tsv"), i=range(config['cores']['quantify_barcodes']), library=LIBRARY, allow_missing=True),
        expand(os.path.join(config['project_dir'], "{library}", "quant_dir", "worker{i}_" + str(config['cores']['quantify_barcodes']) + ".metrics.tsv"), i=range(config['cores']['quantify_barcodes']), library=LIBRARY, allow_missing=True),
        expand(os.path.join(config['project_dir'], "{library}", "quant_dir", "worker{i}_" + str(config['cores']['quantify_barcodes']) + ".ambig.counts.tsv"), i=range(config['cores']['quantify_barcodes']), library=LIBRARY, allow_missing=True),
        expand(os.path.join(config['project_dir'], "{library}", "quant_dir", "worker{i}_" + str(config['cores']['quantify_barcodes']) + ".ambig.partners"), i=range(config['cores']['quantify_barcodes']), library=LIBRARY, allow_missing=True),
        expand(os.path.join(config['project_dir'], "{library}", "{library}.bam"), library=LIBRARY),
        expand(os.path.join(config['project_dir'], "{library}", "{library}.bam.bai"), library=LIBRARY),
        expand(os.path.join(config['project_dir'], "{library}", "{library}.counts.tsv.gz"), library=LIBRARY),
        expand(os.path.join(config['project_dir'], "{library}", "{library}.quant_metrics.tsv.gz"), library=LIBRARY)

rule fastqc_biological_reads:
    input:
        expand(os.path.join("home/adrien.dufour/PROTECT/debug_data/", "{split}_R1.fastq"), split=SPLIT)
    params:
        outdir=os.path.join(config['project_dir'], 'fastqc')
    output:
        [os.path.join(config['project_dir'], 'fastqc', x.replace('.fastq', '_fastqc.html')) for x in FASTQ if 'R1' in x]
    shell:
        "fastqc {input} -o {params.outdir}"

rule filter_reads:
    input:
        fastq=FASTQ,
        yaml=YAMLPATH
    output:
        os.path.join(config['project_dir'], "{library}", "filtered_parts", "{library}_{run}_{split}.fastq"),
        os.path.join(config['project_dir'], "{library}", "filtered_parts", "{library}_{run}_{split}.fastq.counts.pickle"),
        os.path.join(config['project_dir'], "{library}", "filtered_parts", "{library}_{run}_{split}metrics.yaml"),
    conda:
        "/home/adrien.dufour/NeuroDev_ADD/Envs/single_cell.yaml"
    params:
        workers=config['cores']['default']
    log:
        #"logs/{wildcards.run}_{wildcards.library}_{params.worker}_filter.log"
    shell:
        """
        for i in {{0..{params.workers}}}; do
            python indrops.py {input.yaml} filter --runs {RUN} --libraries {LIBRARY} --total-workers {params.workers} --worker-index $i
        done;
        """
        # Resulting workload (a list of run parts), will be split among N --total-workers,
        # where worker with --worker-index i will do steps (i, N+i, 2N+i, ...)

rule abundant_barcodes:
    input:
        expand(os.path.join(config['project_dir'], "{library}", "filtered_parts", "{library}_{run}_{split}.fastq.counts.pickle"), split=SPLIT, library=LIBRARY, run=RUN),
        yaml=YAMLPATH
    output:
        os.path.join(config['project_dir'], "{library}", "abundant_barcodes.pickle"),
        os.path.join(config['project_dir'], "{library}", "{library}.barcode_abundance_by_barcode.png"),
        os.path.join(config['project_dir'], "{library}", "{library}.barcode_abundance.png"),
        os.path.join(config['project_dir'], "{library}", "{library}.filtering_stats.csv")
    conda:
        "/home/adrien.dufour/NeuroDev_ADD/Envs/single_cell.yaml"
    shell:
        """
        python indrops.py {input.yaml} identify_abundant_barcodes --libraries {LIBRARY}
        """

rule sort_reads:
    input:
        expand(os.path.join(config['project_dir'], "{library}", "{library}.filtering_stats.csv"), library=LIBRARY),
        yaml=YAMLPATH
    output:
        os.path.join(config['project_dir'], "{library}", "filtered_parts", "{library}_{run}_{split}.fastq.sorted.fastq.gz"),
        os.path.join(config['project_dir'], "{library}", "filtered_parts", "{library}_{run}_{split}.fastq.sorted.fastq.gz.index.pickle")
    conda:
        "/home/adrien.dufour/NeuroDev_ADD/Envs/single_cell.yaml"
    params:
        workers=config['cores']['default']
    log:
        #"logs/{library}_{split}_sort.log"
    shell:
        """
        for i in {{0..{params.workers}}}; do
            python indrops.py {input.yaml} sort --libraries {LIBRARY} --total-workers {params.workers} --worker-index $i
        done;
        """

rule quantify_barcodes:
    input:
        expand(os.path.join(config['project_dir'], "{library}", "filtered_parts", "{library}_{run}_{split}.fastq.sorted.fastq.gz.index.pickle"), split=SPLIT, run=RUN, library=LIBRARY),
        yaml=YAMLPATH,
    output:
        expand(os.path.join(config['project_dir'], "{library}", "quant_dir", "worker{i}_" + str(config['cores']['quantify_barcodes']) + ".counts.tsv"), i=range(config['cores']['quantify_barcodes']), library=LIBRARY, allow_missing=True),
        expand(os.path.join(config['project_dir'], "{library}", "quant_dir", "worker{i}_" + str(config['cores']['quantify_barcodes']) + ".metrics.tsv"), i=range(config['cores']['quantify_barcodes']), library=LIBRARY, allow_missing=True),
        expand(os.path.join(config['project_dir'], "{library}", "quant_dir", "worker{i}_" + str(config['cores']['quantify_barcodes']) + ".ambig.counts.tsv"), i=range(config['cores']['quantify_barcodes']), library=LIBRARY, allow_missing=True),
        expand(os.path.join(config['project_dir'], "{library}", "quant_dir", "worker{i}_" + str(config['cores']['quantify_barcodes']) + ".ambig.partners"), i=range(config['cores']['quantify_barcodes']), library=LIBRARY, allow_missing=True)
    conda:
        "/home/adrien.dufour/NeuroDev_ADD/Envs/single_cell.yaml"
    params:
        cores=config['cores']['quantify_barcodes'],
        max_idx=config['cores']['quantify_barcodes'] - 1
    shell:
        """
        for i in {{0..{params.cores}}}; do
            python indrops.py {input.yaml} quantify --libraries {LIBRARY} --total-workers {params.cores} --worker-index $i
        done;
        """

rule aggregate_umis:
    input:
        #lambda wildcards: aggregate_input(wildcards),
        expand(os.path.join(config['project_dir'], "{library}", "quant_dir", "worker{i}_" + str(config['cores']['quantify_barcodes']) + ".counts.tsv"), i=range(config['cores']['quantify_barcodes']), library=LIBRARY, allow_missing=True),
        yaml=YAMLPATH
    params:
        workers=config['cores']['default']
    output:
        expand(os.path.join(config['project_dir'], "{library}", "{library}.bam"), library=LIBRARY),
        expand(os.path.join(config['project_dir'], "{library}", "{library}.bam.bai"), library=LIBRARY),
        expand(os.path.join(config['project_dir'], "{library}", "{library}.counts.tsv.gz"), library=LIBRARY),
        expand(os.path.join(config['project_dir'], "{library}", "{library}.quant_metrics.tsv.gz"), library=LIBRARY)
    conda:
        "/home/adrien.dufour/NeuroDev_ADD/Envs/single_cell.yaml"
    shell:
        """
        python indrops.py {input.yaml} aggregate --total-workers {params.workers} --libraries {LIBRARY}
        """