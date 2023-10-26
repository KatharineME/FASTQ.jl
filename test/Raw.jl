using Test: @test

using FASTQ

# ---- #

const TE = FASTQ.TE

const DA = FASTQ._DA

const DAT = joinpath(DA, "Test")

const DAD = joinpath(DAT, "DNA")

const S1 = "Sample1"

# ---- #

const FQ_ = FASTQ.Raw.find(joinpath(DAD, S1))

@test length(FQ_) == 2

@test all([endswith(fq, ".gz") for fq in FQ_])

# ---- #

N_JO = 8

# ---- #

const PAC = FASTQ.Support.trash_remake_directory(joinpath(TE, "CheckRead"))

FASTQ.Raw.check(PAC, FQ_, N_JO)

@test sum(collect(endswith(fi, "fastqc.html") for fi in readdir(PAC))) == 2

# ---- #

const CO = joinpath(TE, "TestConcatenate")

run(`cp -Rf $DAD $CO`)

const S1C = joinpath(CO, S1)

const S2 = "Sample2"

const S2C = replace(S1C, S1 => S2)

run(`cp -Rf $S1C/ $S2C/`)

# ---- #

const FQ2_ = FASTQ.Raw.find(S2C)

FASTQ.Raw.concatenate(CO, FQ2_)

@test round(FASTQ.Support.calculate_size(joinpath(CO, "R1.fastq.gz"))) == 3

# ---- #

const R1 = joinpath(DAT, "DNA", S1, "test_dna_4k.R1.fastq.gz")

const R2 = replace(R1, "R1" => "R2")

const RE_ = (R1, R2)

const PA = joinpath(TE, "CheckRaw")

const SOR1 = joinpath(DAT, "DNA", S2, "test_dna_40k.R1.fastq.gz")

const SOR2 = replace(SOR1, "R1" => "R2")

# ---- #

const GEC, SOC =
    [FASTQ.Support.trash_remake_directory(joinpath(TE, pa)) for pa in ("Germline", "Somatic")]

for (pa, re_, le) in ((GEC, RE_, 2), (SOC, [RE_..., SOR1, SOR2], 4))

    FASTQ.Raw.check(pa, re_, N_JO)

    @test sum(collect(endswith(fi, "fastqc.html") for fi in readdir(pa))) == le

end

# ---- #

const TR = FASTQ.Support.trash_remake_directory(joinpath(TE, "Trim"))

# ---- #

FASTQ.Raw.trim(TR, R1, R2, N_JO)

@test sum([endswith(fi, "fastq.gz") for fi in readdir(TR)]) == 2

@test sum([occursin("fastp", fi) for fi in readdir(TR)]) == 2

# ---- #

const GEN = "GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set"

const GR = joinpath(DA, "GRCh38")

const GEPA = joinpath(GR, GEN)

const GE = joinpath(GEPA, string(GEN, ".fna.gz"))

const ALD = FASTQ.Support.trash_remake_directory(joinpath(TE, "AlignDNA"))

const ME = 8

# ---- #

const BA = FASTQ.Raw.align_dna(ALD, S1, R1, R2, GE, N_JO, ME)

@test round(FASTQ.Support.calculate_size(BA)) > 770

@test sum(occursin("bam", fi) for fi in readdir(ALD)) == 3

# ---- #

const ALCT, ALCG, ALQ, ALS = [
    FASTQ.Support.trash_remake_directory(joinpath(TE, pa)) for pa in (
        "AlignBulkCDNAtoTranscriptome",
        "AlignBulkCDNAtoGenome",
        "AlignandQuantifyBulkCDNAtoGenome",
        "AlignSingleCellCDNAtoGenome",
    )
]

const TRN = "Homo_sapiens.GRCh38.cdna.all"

const TRA = joinpath(GR, TRN, string(TRN, ".fa.gz"))

const CDB = joinpath(DAT, "cDNABulk", S1)

const FR = 51

const SD = 0.05

const R1C = joinpath(CDB, "s1_cdna_4k.R1.fastq.gz")

const R2C = replace(R1C, "R1" => "R2")

# ---- #

FASTQ.Raw.align_bulk_cdna_to_transcriptome(ALCT, R1C, R2C, FR, SD, TRA, N_JO)

@test round(FASTQ.Support.calculate_size(joinpath(ALCT, "abundance.tsv"))) == 6

# ---- #

const BAN = "Aligned.sortedByCoord.out.bam"

const BAT = joinpath(DAT, "DNABAMSTAR", S1)

const BST, BAI = [joinpath(BAT, string(BAN, st)) for st in (".stat", ".bai")]

run(`rm -f $BST $BAI`)

FASTQ.Raw._index_and_stat(BAT, N_JO)

@test all([isfile(fi) for fi in (BST, BAI)])

# ---- #

const ID = FASTQ.Reference.generate_star_genome_file(GE, N_JO; ga = nothing)

# ---- #

const BACG = FASTQ.Raw.align_bulk_cdna_to_genome(ALCG, R1C, R2C, ID, N_JO)

@test round(FASTQ.Support.calculate_size(BACG)) == 791

# ---- #

const BAQ = FASTQ.Raw.align_and_quantify_bulk_cdna_to_genome(ALQ, R1C, R2C, ID, N_JO)

@test round(FASTQ.Support.calculate_size(joinpath(BAQ))) == 791

# ---- #

const DAS = joinpath(DAT, "cDNASingleCell", S1)

const R1S = joinpath(DAS, "400K.R1.fastq.gz")

const R2S = replace(R1S, "R1" => "R2")

# ---- #

FASTQ.Raw.align_single_cell_cdna_to_genome(
    ALS,
    R1S,
    R2S,
    ID,
    N_JO;
    wh = nothing,
    bas = 1,
    bal = 16,
    ums = 17,
    uml = 12,
    rel = 151,
)

@test round(
    FASTQ.Support.calculate_size(joinpath(ALS, "Solo.out", "Gene", "filtered", "matrix.mtx")),
) == 703
