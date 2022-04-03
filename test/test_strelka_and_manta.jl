include("_.jl")

to = joinpath(dirname(@__DIR__), "tool")

Fastq.support.test_strelka_and_manta(to)
