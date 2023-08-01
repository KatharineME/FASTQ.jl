using FASTQ

include("environment.jl")

# ---- #

for pa in (joinpath(TE, "missing_path"), TE)

    BioLab.Path.warn_overwrite(pa)

end

# ---- #

DA = FASTQ.DA

DAT = joinpath(FASTQ.DA, "Test")

DAR = joinpath(DA, "ReferenceGenome", "GRCh38")

DAC = joinpath(DAR, "Chromosome")

TO = joinpath(dirname(@__DIR__), "tool")

# ---- #

n_jo = 8

me = 8

mo = "dna"

rs = true

ex = false

cd = joinpath(DAT, "cDNA")

chs = joinpath(DAC, "chromosome.bed.gz")

chn = joinpath(DAC, "chrn_n.tsv")

sn = joinpath(TO, "snpEff", "snpEff.jar")

va = joinpath(DA, "Ensembl", "GRCh38", "homo_sapiens-grch38-chr1_y.vcf.gz")

ge = joinpath(DAR, "GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz")

# FASTQ.Command.call_variants_on_bulk_cdna(TE, cd, n_jo, ge, mo, ex, chs, chn, me, TO, sn, rs, va)

# ---- #

fr = 51

sd = 0.05

or = "human"

tr = joinpath(DA, "ReferenceTranscriptome", "Homo_sapiens.GRCh38.cdna.all.fa.gz")

mg = joinpath(DA, "Mouse", "mouse_transcript_mouse_gene.tsv")

# FASTQ.Command.measure_gene_expression_of_bulk_cdna(TE, cd, n_jo, tr, fr, sd, or, mg)

# ---- #

# FASTQ.Command.measure_gene_expression_of_single_cell_cdna()

# ---- #

sa = "test"

DAD = joinpath(DA, "Test", "DNA")

r1 = joinpath(DAD, "test_dna_4k.R1.fastq.gz")

r2 = joinpath(DAD, "test_dna_4k.R2.fastq.gz")

# FASTQ.Command.call_variants_on_germline_dna(
#     TE,
#     r1,
#     r2,
#     n_jo,
#     me,
#     sa,
#     ge,
#     chs,
#     chn,
#     sn,
#     mo,
#     ex,
#     TO,
#     rs,
#     va,
# )

# ---- #

# FASTQ.Command.call_variants_on_germline_dna(
#     TE,
#     DAD,
#     n_jo,
#     me,
#     ge,
#     chs,
#     chn,
#     sn,
#     mo,
#     ex,
#     TO,
#     rs,
#     va,
# )

# ---- #

sor1 = joinpath(DAD, "test_dna_40k.R1.fastq.gz")

sor2 = joinpath(DAD, "test_dna_40k.R2.fastq.gz")

# FASTQ.Command.call_variants_on_somatic_dna(
#     TE,
#     r1,
#     r2,
#     sor1,
#     sor2,
#     n_jo,
#     me,
#     sa,
#     ge,
#     chs,
#     chn,
#     sn,
#     mo,
#     ex,
#     TO,
#     rs,
#     va,
# )

# ---- #

rt = joinpath(TO, "rtg-tools-3.11")

vq = "benchmark/apply_germline_dna_to_genome/call_germline_variant/pass.vcf.gz"

vt = "benchmark/HG002_truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"

be = "benchmark/HG002_truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed.gz"

nch = "grch/chromosome/n_chrn.tsv"

# FASTQ.Command.benchmark(TE, ge, rt, n_jo, nch, vq, vt, be)

# ---- #

FASTQ.Command.concatenate_fastq(DAD, "R1")
