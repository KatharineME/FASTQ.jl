include("_.jl")

n_jo, me, mo, _, _, _, sa, _, ou, r1, r2, _, _, ge, _, _, _, _ = Fastq.read_setting(se)

al = joinpath(ou, "align_$mo")

ba = joinpath(al, "$sa.bam")

Fastq.align_dna(al, sa, ba, r1, r2, ge, n_jo, me)
