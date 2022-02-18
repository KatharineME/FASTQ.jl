include("_.jl")

to = joinpath(dirname(@__DIR__), "tool")

Fastq.test_strelka_and_manta(to)
