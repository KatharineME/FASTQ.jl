TE = joinpath(tempdir(), "Fastq.test", "")

if isdir(TE)

    rm(TE; recursive = true)

end

mkdir(TE)

println("Made ", TE)

using Revise
using BenchmarkTools

using Fastq

rm(TE; recursive = true)

println("Removed ", TE)
