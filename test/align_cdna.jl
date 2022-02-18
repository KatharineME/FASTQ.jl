include("_.jl")

n_jo, _, _, _, _, _, sa, _, ou, _, _, _, _, ge, _, _, _, _ = Fastq.read_setting(se)

al = joinpath(ou, "align_cdna")

Fastq.align_cdna(al, sa, cr1, cr2, ge, n_jo)
