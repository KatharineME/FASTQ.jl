function map_mouse_transcript_to_mouse_gene(nu_tr_sa, ma)

    Fastq.support.log()

    mt_mg = CSV.read(ma, DataFrame)

    rename!(mt_mg, "Transcript stable ID" => :id, "Gene name" => :gene)

    dropmissing!(mt_mg, :gene)

    nu_ge_sa = innerjoin(nu_tr_sa, mt_mg, on = :id)

    select!(nu_ge_sa, :gene, Not([:id]))

    nu_ge_sa =
        combine(groupby(nu_ge_sa, :gene), names(nu_ge_sa, Not(:gene)) .=> sum, renamecols = false)

    return nu_ge_sa

end
