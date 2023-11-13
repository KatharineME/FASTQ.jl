using Test: @test

using FASTQ

# ---- #

const TE = FASTQ.TE

const DA = FASTQ._DA

const DAT = joinpath(DA, "Test")

const DAD = joinpath(DAT, "DNA")

# ---- #

const CO = joinpath(TE, "TestConcatenate")

run(`cp -Rf $DAD $CO`)

const S1 = "Sample1"

const S2 = replace(S1, "1" => "2")

const S1C = joinpath(CO, S1)

const S2C = replace(S1C, S1 => S2)

run(`cp -Rf $S1C/ $S2C/`)

# ---- #

FASTQ.Command.concatenate_fastq(CO, read_name_scheme = "R1")

@test sum(endswith(di, "Concatenated") for di in readdir(CO)) == 1

@test round(
    sum([
        FASTQ.Support.calculate_size(fi) for
        fi in readdir(joinpath(CO, string(S2, "Concatenated")), join = true)
    ]),
) == 6

# ---- #

const CON, SEO, SSO, FL =
    [mkdir(joinpath(TE, st)) for st in ("1.Concatenate", "2.Snpeff", "3.Snpsift", "4.Filter")]

const DAV = joinpath(DAT, "VCF", S1)

const VC_ =
    (joinpath(DAV, "Manta", "diploidSV.vcf.gz"), joinpath(DAV, "Strelka", "variants.vcf.gz"))

const DAR = joinpath(DA, "GRCh38")

const GEN = "GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set"

const GEPA = joinpath(DAR, GEN)

const GE = joinpath(GEPA, string(GEN, ".fna.gz"))

const CHN = joinpath(GEPA, "chrn_n.tsv")

const VA = joinpath(DAR, "homo_sapiens-grch38-chr1_y.vcf.gz")

const TO = joinpath(FASTQ.PR, "tool")

const SE = joinpath(TO, "snpEff", "snpEff.jar")

const N_JO = 8

const ME = 8

# ---- #

FASTQ.Command._combine_and_annotate_vcf(CON, SEO, SSO, FL, VC_, GE, CHN, VA, SE, N_JO, ME)

@test round(FASTQ.Support.calculate_size(joinpath(CON, "concat.vcf.gz"))) == 27

@test round(FASTQ.Support.calculate_size(joinpath(SSO, "snpsift.vcf.gz"))) == 150

@test round(FASTQ.Support.calculate_size(joinpath(FL, "pass.vcf.gz"))) == 21

# ---- #

const EX = false

const MO = "dna"

const CHS = joinpath(GEPA, "chromosome.bed.gz")

const AN, SS, SSV = "6.Annotate", "3.Snpsift", "snpsift.vcf.gz"

# ---- #

FASTQ.Command.call_variants_on_germline_dna(TE, DAD, EX, GE, VA, TO, N_JO, ME;)

@test round(
    FASTQ.Support.calculate_size(joinpath(TE, "CallVariantsonGermlineDNA", AN, SS, S1, SSV)),
) == 21

# ---- #

const R1 = joinpath(DAD, S1, "test_dna_4k.R1.fastq.gz")

const R2 = replace(R1, "R1" => "R2")

const SOR1 = joinpath(DAD, S2, "test_dna_40k.R1.fastq.gz")

const SOR2 = replace(SOR1, "R1" => "R2")

# ---- #

@test FASTQ.Command.call_variants_on_somatic_dna(
    TE,
    R1,
    R2,
    SOR1,
    SOR2,
    EX,
    S1,
    GE,
    VA,
    TO,
    N_JO,
    ME,
) == nothing

@test round(FASTQ.Support.calculate_size(joinpath(TE, "CallVariantsonSomaticDNA", AN, SS, SSV))) ==
      19

# ---- #

const CD = joinpath(DAT, "cDNABulk")

# ---- #

@test FASTQ.Command.call_variants_on_bulk_cdna(TE, CD, EX, GE, VA, TO, N_JO, ME) == nothing

@test round(
    FASTQ.Support.calculate_size(
        joinpath(TE, "CallVariantsonBulkCDNA", "4.Annotate", SS, S1, SSV),
    ),
) == 26

# ---- #

const OR = "human"

const FR = 51

const SD = 0.05

const TRN = "Homo_sapiens.GRCh38.cdna.all"

const TR = joinpath(DAR, TRN, string(TRN, ".fa.gz"))

const MG = joinpath(DA, "GRCm38", "mouse_transcript_mouse_gene.tsv")

const MGE = "MeasureGeneExpressionofBulkCDNA"

# ---- #

@test FASTQ.Command.measure_gene_expression_of_bulk_cdna(
    TE,
    CD,
    OR,
    FR,
    SD,
    TR,
    N_JO;
    method = "align_to_transcriptome",
    mouse_transcript_to_mouse_gene = MG,
) == nothing

@test round(
    FASTQ.Support.calculate_size(
        joinpath(TE, MGE, "2.AlignBulkCDNAtoTranscriptome", S1, "abundance.tsv"),
    ),
) == 6

# ---- #

@test FASTQ.Command.measure_gene_expression_of_bulk_cdna(
    TE,
    CD,
    OR,
    FR,
    SD,
    GE,
    N_JO;
    method = "align_to_genome",
    mouse_transcript_to_mouse_gene = MG,
) == nothing

@test round(
    FASTQ.Support.calculate_size(
        joinpath(
            TE,
            MGE,
            "2.AlignandQuantifyBulkCDNAtoGenome",
            S1,
            "Aligned.sortedByCoord.out.bam",
        ),
    ),
) == 791

# ---- #

const SI = joinpath(DAT, "cDNASingleCell")

# ---- #

@test FASTQ.Command.measure_gene_expression_of_single_cell_cdna(
    TE,
    SI,
    GE,
    N_JO;
    gene_annotation = nothing,
    whitelist = nothing,
    barcodestart = 1,
    barcodelength = 16,
    umistart = 17,
    umilength = 12,
    readlength = 151,
) == nothing

@test round(
    FASTQ.Support.calculate_size(
        joinpath(
            TE,
            "MeasureGeneExpressionofSingleCellCDNA",
            "2.AlignSingleCellCDNAtoGenome",
            S1,
            "Solo.out",
            "Gene",
            "filtered",
            "matrix.mtx",
        ),
    ),
) == 703
