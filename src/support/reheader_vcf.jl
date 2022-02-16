function reheader_vcf(sa, pa, n_jo)

    na = string(split(basename(pa), "vcf.gz")[1], "reheader.vcf.gz")

    par = joinpath(dirname(pa), na)

    run(pipeline(`bcftools reheader --threads $n_jo --samples $sa $pa`, "$par"))

    run(`tabix --force $par`)

    return par

end
