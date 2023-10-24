using Test: @test

using FASTQ

# ---- #

const TE = FASTQ.TE

const DA = FASTQ._DA

const S1 = "Sample1"

const BAGE = joinpath(DA, "Test", "DNABAM", S1, string(S1, ".bam"))

const DAR = joinpath(DA, "GRCh38")

const GEN = "GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set"

const GEPA = joinpath(DAR, GEN)

const GE = joinpath(GEPA, string(GEN, ".fna.gz"))

const CHS = joinpath(GEPA, "chromosome.bed.gz")

const TO = joinpath(FASTQ.PR, "tool")

# ---- #

id, voo, vost, vostr, vorfi, voc, vogefi, vot =
    FASTQ.BAM._run_strelka_manta_docker_container(TE, BAGE, GE, CHS, TO)

@test length(readlines(pipeline(`docker exec --interactive $id bash -c "ls /home/GRCh38/"`))) == 22

# ---- #

const ST = "path"

pa_ = FASTQ.BAM._set_output_path(ST)

@test pa_[4] == joinpath(ST, "concat.vcf.gz")

# ---- #

const MO = "dna"

const EX = false

const CHN = joinpath(GEPA, "chrn_n.tsv")

const VA = joinpath(DAR, "homo_sapiens-grch38-chr1_y.vcf.gz")

const N_JO = 8

const ME = 8

# ---- #

@test FASTQ.BAM._set_strelka_manta_run(N_JO, ME) == "--mode local --jobs $N_JO --memGb $ME --quiet"

# ---- #

vc_, paco = FASTQ.BAM.call_germline_variant(TE, BAGE, MO, EX, GE, CHS, CHN, VA, TO, N_JO, ME)

@test length(vc_) == 2

@test all([isfile(fi) for fi in vc_])

# ---- #

const BASO = replace(BAGE, S1 => "Sample2")

# ---- #

vc_, paco = FASTQ.BAM.call_somatic_variant(TE, BAGE, BASO, mo, CHS, CHN, VA, TO, N_JO, ME)

@test length(vc_) == 3

@test all([isfile(fi) for fi in vc_])
