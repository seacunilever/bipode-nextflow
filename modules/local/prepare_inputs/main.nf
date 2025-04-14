process PREPARE_INPUTS {
    publishDir "${results_dir}/", mode: "copy"

    input:
    path meta_data
    path meta_mapper
    path counts
    path substances_cell_types
    val additional_divider
    val batch_key
    val min_percent_mapped_reads
    val min_num_mapped_reads
    val min_avg_treatment_count
    val specific_filters
    val results_dir

    output:
    path "bifrost_inputs/*.json", emit: prepared_inputs

    script:
    """
    prepare_bifrost_inputs.py \
        --meta-data $meta_data \
        --meta-mapper $meta_mapper \
        --counts $counts \
        --substances-cell-types $substances_cell_types \
        --additional-divider "$additional_divider" \
        --batch-key "$batch_key" \
        --min-percent-mapped-reads $min_percent_mapped_reads \
        --min-num-mapped-reads $min_num_mapped_reads \
        --min-avg-treatment-count $min_avg_treatment_count \
        --specific-filters '$specific_filters' \
        --output-dir bifrost_inputs
    """
}
 