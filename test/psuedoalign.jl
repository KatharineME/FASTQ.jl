include("_.jl")

n_jo, _, _, _, fr, sd, _, _, ou, _, _, _, _, _, tr, _, _, _ = Fastq.read_setting(se)

Fastq.psuedoalign(tr, n_jo, ou, cr1, cr2, fr, sd)
