module Abundance

using ..FASTQ

function map_mouse_transcript_to_mouse_gene(nu_tr_sa, ma)

    FASTQ.support.log_sub_level_function()

    mt_mg_ = CSV.read(ma, DataFrame)

    rename!(mt_mg_, "Transcript stable ID" => :id, "Gene name" => :gene)

    dropmissing!(mt_mg_, :gene)

    nu_ge_sa = innerjoin(nu_tr_sa, mt_mg_; on = :id)

    select!(nu_ge_sa, :gene, Not([:id]))

    nu_ge_sa = combine(
        groupby(nu_ge_sa, :gene),
        names(nu_ge_sa, Not(:gene)) .=> sum;
        renamecols = false,
    )

end

function make_gene_by_sample(pap, pou, or, ma)

    FASTQ.Support.log_sub_level_function()

    tpm__ = Dict()

    for (ro, di_, fi_) in walkdir(pap)

        for fi in fi_

            if endswith(fi, ".tsv")

                paf = joinpath(ro, fi)

                ab = CSV.read(paf, DataFrame)

                na = basename(splitdir(paf)[1])

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

    id_ = []

    for id in nu_tr_sa[!, :transcript_id]

        push!(id_, splitext(id)[1])

    end

    nu_tr_sa[!, :id] = id_

    select!(nu_tr_sa, :id, Not(:transcript_id))

    if or == "mouse"

        df = map_mouse_transcript_to_mouse_gene(nu_tr_sa, ma)

        co = :gene

    else

        df = nu_tr_sa

        co = :id

    end

    na_, ma_ = Nucleus.Gene.rename(df[!, co])

    insertcols!(df, 1, :human_gene => na_, :membership => ma_)

    df = df[in.(df.membership, Ref([0, 1])), :]

    select!(df, Not([:membership, co]))

    nu_ge_sa = combine(
        groupby(df, :human_gene),
        names(df, Not(:human_gene)) .=> sum;
        renamecols = false,
    )

    CSV.write(joinpath(pou, "gene_x_sample.tsv"), nu_ge_sa)

end

end
