include("_.jl")

n_jo, me, mo, ta, _, _, sa, to, ou, _, _, _, _, ge, _, chs, chn, sn = Fastq.read_setting(se)

ba = joinpath(ou, "align_dna/$sa.markdup.bam")

pao = joinpath(ou, "call_germline_variant")

Fastq.call_germline_variant(mo, ta, ba, ge, chs, chn, pao, n_jo, me, to, sn)
