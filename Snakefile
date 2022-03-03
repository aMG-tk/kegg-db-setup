import glob
import pandas as pd
from snakemake.utils import validate, min_version


def get_md5sum(x):
    import hashlib

    return hashlib.md5(x.encode("utf-8")).hexdigest()[:10]


"""
Author: A. Fernandez-Guerra
Affiliation: Lundbeck Foundation GeoGenetics Centre
Aim: Assembly and initial analysis of ventilation samples
Run: snakemake   -s Snakefile
"""
#
##### set minimum snakemake version #####
min_version("5.20.1")

# configfile: "config/config.yaml"
# report: "report/workflow.rst"

# This should be placed in the Snakefile.

"""
Working directory
"""


workdir: config["wdir"]


# message("The current working directory is " + WDIR)

"""
 The list of samples to be processed
"""

# files = pd.read_csv(config['sample_file'], sep = '\t')

# #SAMPLES, = glob_wildcards(config['wdir'] + "/data/" + config['exp'] + "/{smp}.fq.gz")
# SAMPLES = list(files.label.unique())
# FILES = list(files.file.unique())

DUPS = ["rmdup"]

sample_table_read = pd.read_table(
    config["sample_file_read"], sep="\t", lineterminator="\n"
)
sample_table_read = sample_table_read.drop_duplicates(
    subset="label", keep="first", inplace=False
)
sample_table_read = sample_table_read.dropna()
sample_table_read["label_md5"] = sample_table_read.apply(
    lambda row: get_md5sum(row.label), axis=1
)
sample_table_read.set_index("label_md5", inplace=True)
sample_label_dict_read = sample_table_read.to_dict()["label"]
sample_label_read = sample_table_read.index.values


rule all:
    input:
        done_initial_stats=expand(
            config["rdir"] + "/stats/{smp}.stats-initial.txt", smp=sample_label_read
        ),
        done_read_extension=expand(
            config["rdir"] + "/read-extension/{smp}.extended.fastq.gz",
            smp=sample_label_read,
        ),
        done_stats_extension=expand(
            config["rdir"] + "/stats/{smp}.stats-extension.txt", smp=sample_label_read
        ),
        done_stats_derep=expand(
            config["rdir"] + "/stats/{smp}.stats-derep.txt", smp=sample_label_read
        ),
        done_stats_derep_summary=(
            config["rdir"] + "/stats/all.stats-derep-summary.tsv.gz"
        ),
        done_stats_initial_summary=(
            config["rdir"] + "/stats/all.stats-initial-summary.tsv.gz"
        ),
        done_stats_extension_summary=(
            config["rdir"] + "/stats/all.stats-extension-summary.tsv.gz"
        ),
        done_read_derep=expand(
            config["rdir"] + "/read-derep/{smp}.fa.gz", smp=sample_label_read
        ),
        done_kraken2_coarse_out=expand(
            config["rdir"] + "/kraken2-profiling-coarse/{smp}.kraken2-coarse.out.gz",
            smp=sample_label_read,
        ),
        done_kraken2_coarse_ids=expand(
            config["rdir"] + "/kraken2-profiling-coarse/{smp}.kraken2-coarse.ids",
            smp=sample_label_read,
        ),
        done_kraken2_coarse_tax=expand(
            config["rdir"]
            + "/kraken2-profiling-coarse/{smp}.kraken2-coarse.tax-summary.tsv.gz",
            smp=sample_label_read,
        ),
        done_kraken2_hires_out=expand(
            config["rdir"] + "/kraken2-profiling-hires/{smp}.kraken2-hires.out.gz",
            smp=sample_label_read,
        ),
        done_kraken2_hires_report=expand(
            config["rdir"]
            + "/kraken2-profiling-hires/{smp}.kraken2-hires-filt.report.gz",
            smp=sample_label_read,
        ),
        done_k2_coarse_summary=(
            config["rdir"] + "/kraken2-profiling-coarse/kraken2-coarse-summary.tsv.gz"
        ),
        done_k2_hires_summary=expand(
            config["rdir"] + "/kraken2-profiling-hires/kraken2-hires-summary.tsv.gz",
            smp=sample_label_read,
        ),
        done_read_ids=expand(
            config["rdir"] + "/kraken2-profiling-hires/{smp}.kraken2-hires.ids",
            smp=sample_label_read,
        ),
        done_k_hires_out_tax=expand(
            config["rdir"]
            + "/kraken2-profiling-hires/{smp}.kraken2-hires.tax-summary.tsv.gz",
            smp=sample_label_read,
        ),
        done_bac_ids=expand(
            config["rdir"] + "/kraken2-profiling-hires/{smp}_Bacteria.txt",
            smp=sample_label_read,
        ),
        done_arc_ids=expand(
            config["rdir"] + "/kraken2-profiling-hires/{smp}_Archaea.txt",
            smp=sample_label_read,
        ),
        done_vir_ids=expand(
            config["rdir"] + "/kraken2-profiling-hires/{smp}_Viruses.txt",
            smp=sample_label_read,
        ),
        dome_noneuk_fastq=expand(
            config["rdir"] + "/read-noneuk/{smp}.read-noneuk.fq.gz",
            smp=sample_label_read,
        ),
        dome_euk_fastq=expand(
            config["rdir"] + "/read-euk/{smp}.read-euk.fq.gz",
            smp=sample_label_read,
        ),
        done_out_none=expand(
            config["rdir"] + "/woltka-profiling/{smp}.woltka-none.tsv.gz",
            smp=sample_label_read,
        ),
        done_out_free=expand(
            config["rdir"] + "/woltka-profiling/{smp}.woltka-free.tsv.gz",
            smp=sample_label_read,
        ),
        done_out_species=expand(
            config["rdir"] + "/woltka-profiling/{smp}.woltka-species.tsv.gz",
            smp=sample_label_read,
        ),
        done_out_genus=expand(
            config["rdir"] + "/woltka-profiling/{smp}.woltka-genus.tsv.gz",
            smp=sample_label_read,
        ),
        done_out_family=expand(
            config["rdir"] + "/woltka-profiling/{smp}.woltka-family.tsv.gz",
            smp=sample_label_read,
        ),
        done_out_order=expand(
            config["rdir"] + "/woltka-profiling/{smp}.woltka-order.tsv.gz",
            smp=sample_label_read,
        ),
        done_out_class=expand(
            config["rdir"] + "/woltka-profiling/{smp}.woltka-class.tsv.gz",
            smp=sample_label_read,
        ),
        done_out_phylum=expand(
            config["rdir"] + "/woltka-profiling/{smp}.woltka-phylum.tsv.gz",
            smp=sample_label_read,
        ),
        done_bowtie2_bam_dedup=expand(
            config["rdir"] + "/woltka-profiling/{smp}.woltka.dedup.bam",
            smp=sample_label_read,
        ),
        done_bowtie2_bam_dedup_filt=expand(
            config["rdir"] + "/woltka-profiling/{smp}.woltka.dedup.filtered.bam",
            smp=sample_label_read,
        ),
        done_bowtie2_bam_dedup_metrics=expand(
            config["rdir"] + "/woltka-profiling/{smp}.woltka.dedup.metrics",
            smp=sample_label_read,
        ),
        done_woltka_stats=expand(
            config["rdir"] + "/woltka-profiling/{smp}.woltka.dedup_stats.tsv.gz",
            smp=sample_label_read,
        ),
        woltka_stats_filt=expand(
            config["rdir"]
            + "/woltka-profiling/{smp}.woltka.dedup_stats-filtered.tsv.gz",
            smp=sample_label_read,
        ),
        done_woltka_summary_none=(
            config["rdir"] + "/woltka-profiling/woltka-summary-none.tsv.gz"
        ),
        done_woltka_summary_free=(
            config["rdir"] + "/woltka-profiling/woltka-summary-free.tsv.gz"
        ),
        done_woltka_summary_species=(
            config["rdir"] + "/woltka-profiling/woltka-summary-species.tsv.gz"
        ),
        done_woltka_summary_genus=(
            config["rdir"] + "/woltka-profiling/woltka-summary-genus.tsv.gz"
        ),
        done_woltka_summary_family=(
            config["rdir"] + "/woltka-profiling/woltka-summary-family.tsv.gz"
        ),
        done_woltka_summary_order=(
            config["rdir"] + "/woltka-profiling/woltka-summary-order.tsv.gz"
        ),
        done_woltka_summary_class=(
            config["rdir"] + "/woltka-profiling/woltka-summary-class.tsv.gz"
        ),
        done_woltka_summary_phylum=(
            config["rdir"] + "/woltka-profiling/woltka-summary-phylum.tsv.gz"
        ),
        done_woltka_summary_map_stats=(
            config["rdir"] + "/woltka-profiling/woltka-summary-mapping.tsv.gz"
        ),
        done_woltka_stats_summary=(
            config["rdir"] + "/woltka-profiling/woltka-mapping.summary.tsv.gz"
        ),
        done_woltka_stats_filtered_summary=(
            config["rdir"] + "/woltka-profiling/woltka-mapping-filtered.summary.tsv.gz"
        ),
        done_mdmg_out_dir=expand(
            config["rdir"] + "/woltka-mdmg/{smp}.woltka-mdmg.weight-{weight}.tar.gz",
            smp=sample_label_read,
            weight=config["mdmg_weight"],
        ),
        done_mdmg_csv=expand(
            config["rdir"] + "/woltka-mdmg/{smp}.woltka-mdmg.weight-{weight}.csv.gz",
            smp=sample_label_read,
            weight=config["mdmg_weight"],
        ),
        done_misincorporation=expand(
            config["rdir"]
            + "/woltka-mdmg/{smp}.woltka-mdmg.misincorporation.weight-{weight}.txt.gz",
            smp=sample_label_read,
            weight=config["mdmg_weight"],
        ),
        done_mdmg_summary=expand(
            config["rdir"] + "/woltka-mdmg/woltka-mdmg.summary.{weight}.tsv.gz",
            weight=config["mdmg_weight"],
        ),


"""
##### load rules #####
"""


include: "rules/stats-initial.smk"
include: "rules/stats-summary-reads.smk"
include: "rules/read-rename.smk"
include: "rules/read-derep.smk"
include: "rules/read-extension.smk"
include: "rules/kraken2-profiling-coarse.smk"
include: "rules/kraken2-profiling-coarse-summary.smk"
include: "rules/kraken2-profiling-get-euk-reads.smk"
include: "rules/kraken2-profiling-hires.smk"
include: "rules/kraken2-profiling-hires-summary.smk"
include: "rules/woltka-profiling.smk"
include: "rules/woltka-profiling-summary.smk"
include: "rules/woltka-metaDMG.smk"
include: "rules/woltka-metaDMG-summary.smk"
