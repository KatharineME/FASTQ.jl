using FASTQ

include("environment.jl")

# ---- #

for pa in (joinpath(TE, "missing_path"), TE)

    BioLab.Path.warn_overwrite(pa)

end

# ---- #

DA = joinpath(FASTQ.DA, "Test")

@test readdir(DA) ==
      ["c2.all.v7.1.symbols.gmt", "gene_x_statistic_x_number.tsv", "h.all.v7.1.symbols.gmt"]

# ---- #

sc_ = [-2.0, -1, 0, 0, 1, 2]

n = length(sc_)

# ---- #
output_directory = TE
read1 =
    read2 =
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
        cdna_read_directory =
            reference_transcriptome = fragment_length,
            fragment_length_standard_deviation,
            organism,
            mouse_transcript_to_mouse_gene,
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
