include("find.jl")

fe_va = Fastq.command.read_setting(se)

Fastq.fastq.check_read(
    re_,
    joinpath(fe_va["output_directory"], "check_read"),
    fe_va["number_of_jobs"],
)
