def get_kos(ko, ko_list, file):
    kos = ko_list.get_group(ko)
    kos = "\n".join(kos["Gene"])
    return kos


rule kegg_genes_cluster:
    input:
        subdb=config["rdir"] + "/kegg-genes-subdb/{ko}/{ko}-db",
        subdb_index=config["rdir"] + "/kegg-genes-subdb/{ko}/{ko}-db.index",
    output:
        fasta=config["rdir"]
        + "/kegg-genes-cluster/{ko}/{ko}-i{identity}-c{coverage}.fasta.gz",
    threads: config["mmseq_cluster_threads"]
    params:
        mmseqs_bin=config["mmseqs_bin"],
        min_seq_id=config["mmseqs_cluster_min_seq_id"],
        coverage=config["mmseqs_cluster_coverage"],
        coverage_mode=config["mmseqs_cluster_coverage_mode"],
        cluster_mode=config["mmseqs_cluster_mode"],
        tmp=config["rdir"] + "/kegg-genes-cluster/{ko}/tmp",
        fasta=config["rdir"]
        + "/kegg-genes-cluster/{ko}/{ko}-i{identity}-c{coverage}.fasta",
        clu_rep=config["rdir"]
        + "/kegg-genes-cluster/{ko}/{ko}-i{identity}-c{coverage}_clu_rep",
        cludb=config["rdir"]
        + "/kegg-genes-cluster/{ko}/{ko}-i{identity}-c{coverage}_clu",
        rdir=config["rdir"] + "/kegg-genes-cluster/{ko}",
        wdir=config["wdir"],
    conda:
        "../envs/concat-db.yaml"
    log:
        config["rdir"] + "/logs/kegg-genes-cluster/{ko}-i{identity}-c{coverage}.log",
    benchmark:
        (
            config["rdir"]
            + "/benchmarks/kegg-genes-cluster/{ko}-i{identity}-c{coverage}.bmk"
        )
    message:
        """--- Cluster KEGG genes database by KO ---"""
    shell:
        """
        echo "--- Cluster KEGG genes ---"
        if [ -d {params.rdir} ]; then
            rm -rf {params.rdir}
            mkdir -p {params.rdir}
        else
            mkdir -p {params.rdir}
        fi

        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        # if index has size larger than zero, then cluster KEGG genes
        N=$(wc -l {input.subdb_index} | awk '{{print $1}}')

        if [[ ${{N}} -gt 0 ]]; then
            {params.mmseqs_bin} cluster {input.subdb} {params.cludb} {params.tmp} \
                -c {params.coverage} \
                --min-seq-id {params.min_seq_id} \
                --cov-mode {params.coverage_mode} \
                --cluster-mode {params.cluster_mode} \
                --threads {threads} 
            {params.mmseqs_bin} result2repseq {input.subdb} {params.cludb} {params.clu_rep} --threads {threads}
            {params.mmseqs_bin} result2flat {input.subdb} {input.subdb} {params.clu_rep} {params.fasta} --use-fasta-header 1
            gzip {params.fasta}
            rm -rf {params.tmp} {params.clu_rep}*
        else
            echo "No KEGG genes to cluster"
            touch {output.fasta}
        fi
        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
