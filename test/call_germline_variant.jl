include("_.jl")

fe_va = Fastq.command.read_setting(se)

ou, sa = fe_va["output_directory"], fe_va["sample"]

ba = joinpath(ou, "align_dna/$sa.markdup.bam")

pao = joinpath(ou, "call_germline_variant")

Fastq.bam.call_germline_variant(
    fe_va["molecule"],
    fe_va["exome"],
    ba,
    fe_va["reference_genome"],
    fe_va["chromosome_position"],
    fe_va["chromosome_name"],
    pao,
    fe_va["number_of_jobs"],
    fe_va["memory"],
    fe_va["tool_directory"],
    fe_va["snpeff"],
)
