module abundance

using CSV
using DataFrames
using BioLab
using ..Fastq

include("make_gene_by_sample.jl")

include("map_mouse_transcript_to_mouse_gene.jl")

end
