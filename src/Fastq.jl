module FASTQ

using Nucleus

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

    if isdir(TE)

        rm(TE; recursive = true)

    end

    mkdir(TE)

end

end
