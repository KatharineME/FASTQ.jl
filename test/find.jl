include("_.jl")

di = joinpath(@__DIR__, "data/dna")

re_ = Fastq.fastq.find(di)
