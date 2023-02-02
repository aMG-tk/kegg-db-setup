def get_kos(ko, ko_list, file):
    kos = ko_list.get_group(ko)
    kos = "\n".join(kos["Gene"])
    return kos


rule kegg_genes_subdb:
    input:
        mmseqs_db=config["rdir"] + "/concat-db/kegg-genes-concat-db",
    output:
        gene_ids=config["rdir"] + "/kegg-genes-subdb/{ko}/{ko}-gene_ids.tsv",
        subdb=config["rdir"] + "/kegg-genes-subdb/{ko}/{ko}-db",
    threads: 2
    params:
        mmseqs_bin=config["mmseqs_bin"],
        kos=lambda wildcards: get_kos(
            wildcards.ko,
            ko_list,
            file=config["rdir"] + "/kegg-genes-subdb/{ko}/{ko}-gene_ids.tsv",
        ),
        rdir=config["rdir"] + "/kegg-genes-subdb/{ko}",
        wdir=config["wdir"],
    conda:
        "../envs/concat-db.yaml"
    log:
        config["rdir"] + "/logs/kegg-genes-subdb/{ko}-subdb.log",
    benchmark:
        config["rdir"] + "/benchmarks/kegg-genes-subdb/{ko}-subdb.bmk"
    message:
        """--- Concatenate KEGG genes databases ---"""
    shell:
        """
        # remove if folder exists
        if [ -d {params.rdir} ]; then
            rm -rf {params.rdir}
            mkdir -p {params.rdir}
        else
            mkdir -p {params.rdir}
        fi
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}
        echo "{params.kos}" > {output.gene_ids}
        {params.mmseqs_bin} createsubdb {output.gene_ids} {input.mmseqs_db} {output.subdb} --id-mode 1 --subdb-mode 1
        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
