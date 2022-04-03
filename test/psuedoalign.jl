include("_.jl")

fe_va = Fastq.command.read_setting(se)

Fastq.fastq.psuedoalign(
    fe_va["reference_transcriptome"],
    fe_va["number_of_jobs"],
    fe_va["output_directory"],
    fe_va["cdna_read1"],
    fe_va["cdna_read2"],
    fe_va["fragment_length"],
    fe_va["fragment_length_standard_deviation"],
)
