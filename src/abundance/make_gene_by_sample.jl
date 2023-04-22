function make_gene_by_sample(pap, pou, or, ma)

    Fastq.support.log()

    # Make transcript by sample from adbundance.tsvs

    tpm__ = Dict()

    for (ro, di_, fi_) in walkdir(pap)

        for fi in fi_
            if endswith(fi, ".tsv")

                paf = joinpath(ro, fi)

                ab = CSV.read(paf, DataFrame;)

                na = splitdir(splitdir(paf)[1])[2]

                tpm__[na] = ab[!, "tpm"]

                if !issubset(["transcript_id"], keys(tpm__))

                    tpm__["transcript_id"] = ab[!, "target_id"]

                end

            end

        end

    end

    nu_tr_sa = DataFrame(tpm__)

    select!(nu_tr_sa, :transcript_id, Not([:transcript_id]))

    CSV.write(joinpath(pou, "transcript_x_sample.tsv"), nu_tr_sa)


    # Remove transcript version number

    id_ = []

    for id in nu_tr_sa[!, :transcript_id]

        push!(id_, splitext(id)[1])

    end

    nu_tr_sa[!, :id] = id_

    select!(nu_tr_sa, :id, Not(:transcript_id))


    # If mouse, convert mouse transcript to mouse gene

    if or == "mouse"

        df = map_mouse_transcript_to_mouse_gene(nu_tr_sa, ma)

        co = :gene

    else

        df = nu_tr_sa

        co = :id

    end


    # Convert to human gene

    na_, ma_ = Biolab.Gene.rename(df[!, co])

    insertcols!(df, 1, :human_gene => na_, :membership => ma_)

    df = df[in.(df.membership, Ref([0, 1])), :]

    select!(df, Not([:membership, co]))

    nu_ge_sa =
        combine(groupby(df, :human_gene), names(df, Not(:human_gene)) .=> sum, renamecols = false)

    CSV.write(joinpath(pou, "gene_x_sample.tsv"), nu_ge_sa)

end
