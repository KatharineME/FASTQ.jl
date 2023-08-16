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

readdir(joinpath(homedir(), ".Trash", basename(DI))) == ["Applejack.txt", "RainbowDash.txt"]

# ---- #
# ---- #
# ---- #
