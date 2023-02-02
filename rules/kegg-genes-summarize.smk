from itertools import chain
from functools import reduce
from pathlib import Path
import shutil


def is_non_zero_file(fpath):
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0


def get_fasta_files(wildcards, identity, coverage):
    files = expand(
        config["rdir"]
        + "/kegg-genes-cluster/{ko}/{ko}-i"
        + identity
        + "-c"
        + coverage
        + ".fasta.gz",
        ko=kos,
    )
    return files


rule kegg_genes_summarize:
    input:
        fasta_files=lambda wc: get_fasta_files(
            wc, identity=wc.identity, coverage=wc.coverage
        ),
    output:
        file=config["rdir"]
        + "/kegg-genes-summary/kegg_genes-i{identity}-c{coverage}.fasta.gz",
    threads: 1
    log:
        config["rdir"]
        + "/logs/kegg-genes-summary/kegg_genes-i{identity}-c{coverage}.log",
    benchmark:
        (
            config["rdir"]
            + "/benchmarks/kegg-genes-summary/kegg_genes-i{identity}-c{coverage}.bmk"
        )
    message:
        """--- Concatenate clustered KEGG fasta files ---"""
    run:
        with open(output.file, "wb") as wfp:
            for fn in input.fasta_files:
                if is_non_zero_file(fn):
                    with open(fn, "rb") as rfp:
                        shutil.copyfileobj(rfp, wfp)
