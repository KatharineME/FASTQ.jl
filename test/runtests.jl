using Test

using FASTQ

FASTQ.Command.call_variants_on_bulk_cdna(
    output_directory,
    cdna_read_directory,
    number_of_jobs,
    reference_genome,
    molecule,
    exome,
    chromosome_position,
    chromosome_name,
    memory,
    tool_directory,
    snpeff,
    annotate_with_rsid,
    variant_database,
)

FASTQ.Command.measure_gene_expression_of_bulk_cdna(
    output_directory,
    cdna_read_directory,
    number_of_jobs,
    reference_transcriptome,
    fragment_length,
    fragment_length_standard_deviation,
    organism,
    mouse_transcript_to_mouse_gene,
)

FASTQ.Command.measure_gene_expression_of_single_cell_cdna()

FASTQ.Command.call_variants_on_germline_dna(
    output_directory,
    read1,
    read2,
    number_of_jobs,
    memory,
    sample,
    reference_genome,
    chromosome_position,
    chromosome_name,
    snpeff,
    molecule,
    exome,
    tool_directory,
    annotate_with_rsid,
    variant_database,
)

FASTQ.Command.call_variants_on_somatic_dna(
    output_directory,
    read1,
    read2,
    somatic_read1,
    somatic_read2,
    number_of_jobs,
    memory,
    sample,
    reference_genome,
    chromosome_position,
    chromosome_name,
    snpeff,
    molecule,
    exome,
    tool_directory,
    annotate_with_rsid,
    variant_database,
)

FASTQ.Command.benchmark(
    output_directory,
    reference_genome,
    rtg_tools,
    number_of_jobs,
    name_chromosome,
    query_vcf,
    truth_vcf,
    confident_regions_bed,
)

FASTQ.Command.concatenate_fastq(dna_read_directory, read_name_scheme)

@info "Tests passed."
