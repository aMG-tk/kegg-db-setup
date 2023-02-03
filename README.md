# A Snakemake workflow to create a dereplicated KEGG database

A simple workflow that creates a dereplicated KEGG GENES database by KO. The workflow needs a KEGG GENES, only available by [subscription](https://www.pathway.jp/en/academic.html) at the moment.

The workflow performs the following steps:
1. Map the KEGG GENES to KEGG Orthology (KO) and create a global MMseqs2 database for the KEGG GENES.
2. Create a MMseqs2 subdatabase for each KO containing the KEGG GENES.
3. Dereplicate the KEGG GENES for each KO using the `cluster` module of MMseqs2.
4. Get representative sequences for each KO using the `result2repseq` module of MMseqs2 and create a fastA file for each KO.
5. Combine all fastA files into a single fastA file.

One can run the workflow using the following command:

```bash
snakemake --snakefile /vol/cloud/geogenetics/repos/kegg-db-setup/Snakefile -d ./ \
    --configfile config/config.yaml --use-conda -j 100 \
    --conda-frontend mamba --latency-wait 60 \
    --cluster-config config/cluster.yaml \
    --cluster "sbatch --export=ALL -t {cluster.time} --ntasks-per-node {cluster.ntasks_per_node} --nodes {cluster.nodes} --cpus-per-task {cluster.cpus_per_task} --partition {cluster.partition} --job-name {rulename}.{jobid} --output=$(pwd)/slurm-%j.out" 
```

> Example using SLURM as a job scheduler.

In the [config](./config/config.yaml) file, one would be interested in modifying the following parameters:
```yaml
# Clustering parameters
mmseqs_cluster_min_seq_id: 0.9
mmseqs_cluster_coverage: 0.8
mmseqs_cluster_coverage_mode: 0
mmseqs_cluster_mode: 2
```