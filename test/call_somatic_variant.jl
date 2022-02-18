include("_.jl")

n_jo, me, _, ta, _, _, sa, to, ou, _, _, _, _, ge, _, chs, chn, sn = Fastq.read_setting(se)

ger = joinpath(ou, "align_dna/$sa.markdup.bam")

som = joinpath(ou, "align_cdna/$(sa).Aligned.sortedByCoord.out.bam")

pao = joinpath(ou, "call_somatic_variant")

Fastq.call_somatic_variant(ta, ger, som, ge, chs, chn, pao, n_jo, me, to, sn)
