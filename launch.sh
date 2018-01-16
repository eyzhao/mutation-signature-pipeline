source /projects/ezhao_prj/dependencies/miniconda3/bin/activate dependencies

if [ "$HOSTNAME" == "n104" ]; then
    echo 'Running on numbers cluster'
    snakemake -p \
        --cluster-config "config/cluster.json" \
        --drmaa ' --mem-per-cpu={cluster.mem_mb} {cluster.flags} ' \
        --jobs 500 \
        --latency-wait 60 \
        --max-jobs-per-second 1 \
        --restart-times 5 \
        --rerun-incomplete \
        --keep-going

    mkdir -p logs/cluster_logs
    mv slurm*.out logs/cluster_logs

elif [ "$HOSTNAME" == "login-apollo.hpc.bcgsc.ca" ]; then
    echo "Running on APOLLO cluster"
    mkdir -p logs/sge_logs_error logs/sge_logs_output
    snakemake -p \
        --cluster-config "config/cluster.json" \
        --cluster "qsub -q arc.q -P arc.prj -V -N msig_timing -e logs/sge_logs_error -o logs/sge_logs_output -l h_vmem={cluster.mem_mb}M" \
        --jobs 500 \
        --latency-wait 60 \
        --max-jobs-per-second 1 \
        --restart-times 5 \
        --rerun-incomplete

else
    echo "Running in non-cluster mode"

    snakemake -p \
        --cores 100 \
        --latency-wait 60 \
        --max-jobs-per-second 1 \
        --restart-times 5 \
        --rerun-incomplete
fi
