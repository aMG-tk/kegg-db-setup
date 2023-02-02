import glob
import pandas as pd
import numpy as np
import tarfile
import io
from snakemake.utils import validate, min_version
import datatable as dt

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

ko_list_file = config["sdir"] + "/ko/ko_genes.list"
ko_list = dt.fread(
    ko_list_file,
    sep="\t",
    header=False,
    nthreads=config["seqkit_threads"],
    columns=["KO", "Gene"],
).to_pandas()
#ko_list = pd.read_csv(ko_list_file, sep="\t", names=["KO", "Gene"])
#ko_list.replace("ko:", "")
ko_list["KO"] = ko_list["KO"].str.replace("ko:", "")
ko_list = ko_list.groupby("KO")
kos = [ko for ko, _ in ko_list]

localrules: kegg_genes_summarize

rule all:
    input:
        done_concat_db=config["rdir"] + "/concat-db/kegg-genes-concat-db.aa.fa.gz",
        done_mmseqs_db=config["rdir"] + "/concat-db/kegg-genes-concat-db",
        # done_cluster_fasta=expand(
        #     config["rdir"] + "/kegg-genes-subdb/{ko}/{ko}.aa.fa.gz", ko=kos
        # ),
        done_gene_ids=expand(
            config["rdir"] + "/kegg-genes-subdb/{ko}/{ko}-gene_ids.tsv", ko=kos
        ),
        done_subdb=expand(config["rdir"] + "/kegg-genes-subdb/{ko}/{ko}-db", ko=kos),
        done_fasta=expand(
            config["rdir"]
            + "/kegg-genes-cluster/{ko}/{ko}-i{identity}-c{coverage}.fasta.gz",
            ko=kos,
            identity=config["mmseqs_cluster_min_seq_id"],
            coverage=config["mmseqs_cluster_coverage"],
        ),
        done_summary_files=expand(
            config["rdir"]
            + "/kegg-genes-summary/kegg_genes-i{identity}-c{coverage}.fasta.gz",
            identity=config["mmseqs_cluster_min_seq_id"],
            coverage=config["mmseqs_cluster_coverage"],
        ),


"""
##### load rules #####
"""


include: "rules/kegg-genes-concatenate.smk"
include: "rules/kegg-genes-subdb.smk"
include: "rules/kegg-genes-cluster.smk"
include: "rules/kegg-genes-summarize.smk"
