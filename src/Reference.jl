module Reference

using ..FASTQ

function get_chromosome_file_path(re)

    gr = dirname(re)

    chn = joinpath(gr, "chrn_n.tsv")

    chs = joinpath(gr, "chromosome.bed.gz")

    chs, chn

end

function index_genome_file(ge, chs)

    geu = ge

    if endswith(ge, ".gz")

        un = splitext(ge)[1]

        if isfile(un)

            geu = un

        end

    end

    for fi in (ge, geu)

        if !(isfile("$fi.fai") && ispath("$fi.gzi"))

            run(`samtools faidx $fi`)

        end

    end

    if !ispath("$chs.tbi")

        run(`tabix --force $chs`)

    end

    nothing

end

function generate_star_genome_file(ge, n_jo; ga = nothing)

    id = joinpath(dirname(ge), "StarIndex")

    if ga == nothing

        ga = joinpath(
            FASTQ._DA,
            "GRCh38",
            "GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set",
            "GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf",
        )

    end

    if !ispath(joinpath(id, "geneInfo.tab"))

        if ispath(id)

            rm(id, recursive = true)

        end

        mkdir(id)

        ged = splitext(ge)[1]

        if !isfile(ged)

            run(pipeline(`bgzip --decompress --stdout $ge`, ged))

        end

        run(
            `star --runMode genomeGenerate --runThreadN $n_jo --genomeDir $id --genomeFastaFiles $ged --sjdbGTFfile $ga`,
        )

    end

    id

end

end
