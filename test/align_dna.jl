include("_.jl")

fe_va = Fastq.command.read_setting(se)

sa = fe_va["sample"]

al = joinpath(fe_va["output_directory"], "align_$(fe_va["molecule"])")

ba = joinpath(al, "$sa.bam")

Fastq.fastq.align_dna(
    al,
    sa,
    ba,
    fe_va["read1"],
    fe_va["read2"],
    fe_va["reference_genome"],
    fe_va["number_of_jobs"],
    fe_va["memory"],
)
