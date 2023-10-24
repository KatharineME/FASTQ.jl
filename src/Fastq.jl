module FASTQ

using BioLab

const PR = dirname(@__DIR__)

const TE = joinpath(tempdir(), "FASTQ")

const _DA = joinpath(PR, "data")

const _MA = "manta-1.6.0.centos6_x86_64"

const _ST = "strelka-2.9.10.centos6_x86_64"

for jl in readdir(@__DIR__)

    if !startswith(jl, '_') && jl != "FASTQ.jl"

        include(jl)

    end

end

function __init__()

    BioLab.Path.remake_directory(TE)

end

end
