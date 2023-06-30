using FASTQ

include("environment.jl")

# ---- #

for pa in (joinpath(TE, "missing_path"), TE)

    BioLab.Path.warn_overwrite(pa)

end

# ---- #

DA = FASTQ.DA

DAT = joinpath(FASTQ.DA, "Test")

DAC = joinpath(DAT, "ReferenceGenome", "Chromosome")

TO = joinpath(dirname(@__DIR__), "tool")

# ---- #

n_jo = 8

me = 8G

ge = joinpath(
    DA,
    "ReferenceGenome",
    "GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz",
)

chs = joinpath(DAC, "chromosome.bed.gz")

chn = joinpath(DAC, "chrn_n.tsv")

cd = joinpath(DA, "cDNA"),
sn = joinpath(TO, "snpEff", "snpEff.jar")

va = joinpath(DA, "Ensembl", "homo_sapiens-chr1_y.vcf.gz")

# ---- #

FASTQ.Command.call_variants_on_bulk_cdna(
    TE,
    cd,
    n_jo,
    ge,
    "dna",
    false,
    chs,
    chn,
    me,
    TO,
    sn,
    true,
    va,
)

# FASTQ.Command.measure_gene_expression_of_bulk_cdna(
#     TE,
#     cdna_read_directory,
#     number_of_jobs,
#     reference_transcriptome,
#     fragment_length,
#     fragment_length_standard_deviation,
#     organism,
#     mouse_transcript_to_mouse_gene,
# )
# 
# FASTQ.Command.measure_gene_expression_of_single_cell_cdna()
# 
# FASTQ.Command.call_variants_on_germline_dna(
#     output_directory,
#     read1,
#     read2,
#     number_of_jobs,
#     memory,
#     sample,
#     reference_genome,
#     chromosome_position,
#     chromosome_name,
#     snpeff,
#     molecule,
#     exome,
#     tool_directory,
#     annotate_with_rsid,
#     variant_database,
# )
# 
# FASTQ.Command.call_variants_on_somatic_dna(
#     output_directory,
#     read1,
#     read2,
#     somatic_read1,
#     somatic_read2,
#     number_of_jobs,
#     memory,
#     sample,
#     reference_genome,
#     chromosome_position,
#     chromosome_name,
#     snpeff,
#     molecule,
#     exome,
#     tool_directory,
#     annotate_with_rsid,
#     variant_database,
# )
# 
# FASTQ.Command.benchmark(
#     output_directory,
#     reference_genome,
#     rtg_tools,
#     number_of_jobs,
#     name_chromosome,
#     query_vcf,
#     truth_vcf,
#     confident_regions_bed,
# )
# 
# FASTQ.Command.concatenate_fastq(dna_read_directory, read_name_scheme)
