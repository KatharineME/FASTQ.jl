using Test: @test

using FASTQ

# ---- #

const TE = FASTQ.TE

const DA = FASTQ._DA

const VC = FASTQ.Support.trash_remake_directory(joinpath(TE, "VCF"))

const DAT = joinpath(DA, "Test")

const DAV = joinpath(DAT, "VCF", "Sample1")

const VC_ = (
    joinpath(DAV, "Manta", "diploidSV.vcf.gz"),
    joinpath(DAV, "Strelka", "variants.vcf.gz"),
)

const DAR = joinpath(DA, "GRCh38")

const GEN = "GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set"

const GEPA = joinpath(DAR, GEN)

const CHN = joinpath(GEPA, "chrn_n.tsv")

const N_JO = 8

# ---- #

co = FASTQ.VCF.combine_vcf(VC, VC_, CHN, N_JO)

@test round(FASTQ.Support.calculate_size(co)) == 27

# ---- #

const VCS = joinpath(DAT, "VCFSomatic")

const SA = joinpath(VCS, "Sample.txt")

# ---- #

par = FASTQ.VCF.reheader_vcf(joinpath(VCS, "somatic.indels.vcf.gz"), SA, N_JO)

@test split(readchomp(`tabix -H $par`), "\t")[10:11] == ["Germline", "Somatic"]

# ---- #

const GE = joinpath(GEPA, string(GEN, ".fna.gz"))

# ---- #

@test FASTQ.VCF._check_genome_version(GE) == 38

# ---- #

const PASE = FASTQ.Support.trash_remake_directory(joinpath(VC, "Snpeff"))

const SE = joinpath(FASTQ.PR, "tool", "snpEff", "snpEff.jar")

const ME = 8

# ---- #

vcse = FASTQ.VCF.annotate_with_snpeff(PASE, co, GE, SE, N_JO, ME)

@test round(FASTQ.Support.calculate_size(vcse)) == 139

# ---- #

const PASS = FASTQ.Support.trash_remake_directory(joinpath(VC, "Snpsift"))

const VA = joinpath(DAR, "homo_sapiens-grch38-chr1_y.vcf.gz")

# ---- #

vcss = FASTQ.VCF.annotate_with_snpsift(PASS, vcse, VA, SE, N_JO)

@test round(FASTQ.Support.calculate_size(vcss)) == 150

# ---- #

vcpa = FASTQ.VCF.filter_vcf(VC, vcss, N_JO)

@test round(FASTQ.Support.calculate_size(vcpa)) == 21
