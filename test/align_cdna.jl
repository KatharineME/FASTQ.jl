include("_.jl")

fe_va = Fastq.command.read_setting(se)

al = joinpath(fe_va["output_directory"], "align_cdna")

Fastq.fastq.align_cdna(
    al,
    fe_va["sample"],
    fe_va["cdna_read1"],
    fe_va["cdna_read2"],
    fe_va["reference_genome"],
    fe_va["number_of_jobs"],
)
