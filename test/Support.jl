using Test: @test

using FASTQ

# ---- #

DI = joinpath(FASTQ.TE, "Ponies")

mkdir(DI)

FI_ = (("RainbowDash.txt", "rainbow"), ("Applejack.txt", "apple"))

for (fi, tx) in FI_

    write(joinpath(DI, fi), tx)

end

FASTQ.Support.trash_remake_directory(DI)

@test isempty(readdir(DI))

@test readdir(joinpath(homedir(), ".Trash", basename(DI))) == ["Applejack.txt", "RainbowDash.txt"]

# ---- #

DII = joinpath(FASTQ.TE, "index")

mkdir(DII)

GR = joinpath(FASTQ._DA, "ReferenceGenome", "GRCh38")

GE = joinpath(GR, "GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz")

CHS = joinpath(GR, "Chromosome", "chromosome.bed.gz")

cp(GE, joinpath(DII, basename(GE)))

cp(CHS, joinpath(DII, basename(CHS)))

FASTQ.Support.index_genome_files(joinpath(DII, basename(GE)), joinpath(DII, basename(CHS)))

# ---- #
