using Test: @test

using FASTQ

# ---- #

const TE = FASTQ.TE

const DA = FASTQ._DA

const S1 = "Sample1"

const BAGE = joinpath(DA, "Test", "DNABAMMinimap", S1, string(S1, ".bam"))

const DAR = joinpath(DA, "GRCh38")

const GEN = "GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set"

const GEPA = joinpath(DAR, GEN)

const GE = joinpath(GEPA, string(GEN, ".fna.gz"))

const CHS = joinpath(GEPA, "chromosome.bed.gz")

const TO = joinpath(FASTQ.PR, "tool")

# ---- #

id = FASTQ.BAM._run_strelka_manta_docker_container(TE, BAGE, GE, CHS, TO)[1]

@test lastindex(
    readlines(pipeline(`docker exec --interactive $id bash -c "ls /home/$GEN/"`)),
) >= 11

FASTQ.Support.remove_docker_container(id)

# ---- #

pa_ = FASTQ.BAM._set_output_path("path")

@test pa_[3] == joinpath("results", "variants")

# ---- #

const MO = "dna"

const EX = false

const CHN = joinpath(GEPA, "chrn_n.tsv")

const VA = joinpath(DAR, "homo_sapiens-grch38-chr1_y.vcf.gz")

const N_JO = 8

const ME = 8

# ---- #

@test FASTQ.BAM._set_strelka_manta_run(N_JO, ME) ==
      "--mode local --jobs $N_JO --memGb $ME --quiet"

# ---- #

const TEG = mkdir(joinpath(TE, "CallGermlineVariant"))

# ---- #

vcg_ = FASTQ.BAM.call_germline_variant(TEG, BAGE, MO, EX, GE, CHS, TO, N_JO, ME)

@test lastindex(vcg_) == 2

@test all([isfile(fi) for fi in vcg_])

# ---- #

const BASO = replace(BAGE, S1 => "Sample2")

const TES = mkdir(joinpath(TE, "CallSomaticVariant"))

# ---- #

vcs_ = FASTQ.BAM.call_somatic_variant(TES, BAGE, BASO, EX, GE, CHS, TO, N_JO, ME)

@test lastindex(vcs_) == 3

@test all([isfile(fi) for fi in vcs_])
