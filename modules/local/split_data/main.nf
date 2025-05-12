process SPLIT_DATA {
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/33/33499ba0fef01239be35b4b8ffae2f35bc921bd88d622a5f5f5c6ed2edb3eaa0/data' :
        'community.wave.seqera.io/library/python_cmdstanpy_numpy_pandas_pruned:b21b7854a692918a' }"

    input:
    tuple val(meta), path(input_data)

    output:
    tuple val(meta), path("${prefix}_batch*.manifest.csv"), path("${prefix}_batch*.tar.gz"), emit: probe_files

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''

    def batch_size = (args2 =~ /--batch-size\s+(\d+)/) ? (args2 =~ /--batch-size\s+(\d+)/)[0][1].toInteger() : 0
    def batch_mode = (args2 =~ /--batch-mode\s+(batch|all)/) ? (args2 =~ /--batch-mode\s+(batch|all)/)[0][1] : 'all'

    """
    split_data.py --input-file $input_data --analysis-dir .

    sleep 5

    cd Data/
    # Get list of all files
    files=(\$(ls *.pkl))
    total_files=\${#files[@]}

    # Create single tar.gz if batch_mode is 'all' or batch_size is 0
    if [ "$batch_mode" = "all" ]; then
        tar -czf "../${prefix}_batch0.tar.gz" \${files[@]}
    fi

    # Split files into batches for manifest creation
    batch_num=1
    for ((i=0; i<total_files; i+=$batch_size)); do
        batch_files=(\${files[@]:i:$batch_size})
        batch_prefix="../${prefix}_batch\${batch_num}"
        manifest_file="\${batch_prefix}.manifest.csv"

        # Create manifest file for this batch
        for f in "\${batch_files[@]}"; do
            echo "\${f%.pkl}" >> \$manifest_file
        done

        if [ "$batch_mode" = "batch" ]; then
            # Create individual tar file for this batch
            tar -czf "\${batch_prefix}.tar.gz" \${batch_files[@]}
        fi

        ((batch_num++))
    done
    """
}
