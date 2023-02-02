rule kegg_genes_concatenate:
    input:
        prok=config["sdir"] + "/fasta/prokaryotes.pep.gz",
        euk=config["sdir"] + "/fasta/eukaryotes.pep.gz",
        t10000=config["sdir"] + "/fasta/T10000.pep.gz",
        t40000=config["sdir"] + "/fasta/T40000.pep.gz",
    output:
        concat_db=config["rdir"] + "/concat-db/kegg-genes-concat-db.aa.fa.gz",
        mmseqs_db=config["rdir"] + "/concat-db/kegg-genes-concat-db",
    threads: config["seqkit_threads"]
    params:
        seqkit_bin=config["seqkit_bin"],
        mmseqs_bin=config["mmseqs_bin"],
        rdir=config["rdir"] + "/concat-db",
        wdir=config["wdir"],
    conda:
        "../envs/concat-db.yaml"
    log:
        config["rdir"] + "/logs/kegg-genes-concat-db/kegg-genes-concat-db.log",
    benchmark:
        config["rdir"] + "/benchmarks/kegg-genes-concat-db/kegg-genes-concat-db.bmk"
    message:
        """--- Concatenate KEGG genes databases ---"""
    shell:
        """
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        if [ -f {output.concat_db} ]; then
            rm -f {output.concat_db}
        fi

        cat {input.prok} {input.euk} {input.t10000} {input.t40000} | {params.seqkit_bin} replace -j {threads} -p "\s.+" -o {output.concat_db}

        N_concat=$({params.seqkit_bin} stats -T -j {threads} {output.concat_db} | tail -n+2 | cut -f 4)
        N_prok=$({params.seqkit_bin} stats -T -j {threads} {input.prok} | tail -n+2 | cut -f 4)
        N_euk=$({params.seqkit_bin} stats -T -j {threads} {input.euk} | tail -n+2 | cut -f 4)
        N_t10000=$({params.seqkit_bin} stats -T -j {threads} {input.t10000} | tail -n+2 | cut -f 4)
        N_t40000=$({params.seqkit_bin} stats -T -j {threads} {input.t40000} | tail -n+2 | cut -f 4)
        N_total=$(echo "scale=0; $N_prok + $N_euk + $N_t10000 + $N_t40000" | bc)

        if [ ${{N_concat}} -ne ${{N_total}} ]; then
            echo "ERROR: concatenated database has $N_concat sequences, but the total number of sequences is $N_total"
            exit 1
        fi

        {params.mmseqs_bin} createdb {output.concat_db} {output.mmseqs_db}

        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
