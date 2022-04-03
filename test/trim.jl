include("_.jl")

fe_va = Fastq.command.read_setting(se)

Fastq.fastq.trim(
    fe_va["read1"],
    fe_va["read2"],
    joinpath(fe_va["output_directory"], "trim"),
    fe_va["number_of_jobs"],
)
