include("_.jl")

fe_va = Fastq.command.read_setting(se)

ou, sa = fe_va["output_directory"], fe_va["sample"]

geba = joinpath(ou, "align_dna/$sa.markdup.bam")

soba = joinpath(ou, "align_cdna/$sa.Aligned.sortedByCoord.out.bam")

pao = joinpath(ou, "call_somatic_variant")

Fastq.bam.call_somatic_variant(
    fe_va["exome"],
    geba,
    soba,
    fe_va["reference_genome"],
    fe_va["chromosome_position"],
    fe_va["chromosome_name"],
    pao,
    fe_va["number_of_jobs"],
    fe_va["memory"],
    fe_va["tool_directory"],
    fe_va["snpeff"],
)
