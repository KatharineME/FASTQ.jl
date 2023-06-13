module Abundance

using CSV

using DataFrames

using BioLab

using FASTQ

function map_mouse_transcript_to_mouse_gene(nu_tr_sa, ma)

    FASTQ.support.log()

    mt_mg = CSV.read(ma, DataFrame)

    rename!(mt_mg, "Transcript stable ID" => :id, "Gene name" => :gene)

    dropmissing!(mt_mg, :gene)

    nu_ge_sa = innerjoin(nu_tr_sa, mt_mg, on = :id)

    select!(nu_ge_sa, :gene, Not([:id]))

    nu_ge_sa =
        combine(groupby(nu_ge_sa, :gene), names(nu_ge_sa, Not(:gene)) .=> sum, renamecols = false)

    return nu_ge_sa

end


pas = "/Users/kate/Downloads/ferreira_treg/output/7166-MR-100/"

# Use the raw STARsolo output. If false, the filtered output will be used.
ra = true

or = "human"

using BioLab

using CSV

using DataFrames

function make_gene_by_sample(pas, fi, ma)

    if ra == true

        pa = "Solo.out/Gene/raw/"

    else

        pa = "Solo.out/Gene/filtered/"

    end

    # ma = CSV.read(joinpath(pas, pa, "matrix.mtx"), DataFrame; delim='\t')

    # ce = CSV.read(joinpath(pas, pa, "barcodes.tsv"), DataFrame; delim='\t')

    # ge = CSV.read(joinpath(pas, pa, "features.tsv"), DataFrame; delim='\t')

    ma = DataFrame(feature = [1, 5, 3], cell = [1, 1, 2], count = [1, 15, 7])

    ge_ = ["A", "B", "C", "D", "E"]

    ce_ = ["cell1", "cell2", "cell3"]

    df = DataFrame()

    df.gene = ma

    # Build gene by cell

    # Plot expression per gene conveniently 

    # Plot mitochondrial expression conveniently

end


function make_gene_by_sample(pap, pou, or, ma)

    FASTQ.support.log()

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

    na_, ma_ = BioLab.Gene.rename(df[!, co])

    insertcols!(df, 1, :human_gene => na_, :membership => ma_)

    df = df[in.(df.membership, Ref([0, 1])), :]

    select!(df, Not([:membership, co]))

    nu_ge_sa =
        combine(groupby(df, :human_gene), names(df, Not(:human_gene)) .=> sum, renamecols = false)

    CSV.write(joinpath(pou, "gene_x_sample.tsv"), nu_ge_sa)

end

end
