module Raw

using ..FASTQ

function find(pa)

    fi_ = readdir(pa, join = true)

    fq_ = fi_[findall(occursin.(".gz", fi_))]

end

function check(pa, fq_, n_jo)

    FASTQ.Support.log_sub_level_function()

    th = minimum((length(fq_), n_jo))

    run(`fastqc --threads $th --outdir $pa $fq_`)

    #run(`multiqc --outdir $pa $pa`)

    basename(pa), fq_[1], fq_[2]

end

function concatenate(pa, fq_; na = "R1")

    FASTQ.Support.log_sub_level_function()

    fo_ = Vector{String}()

    re_ = Vector{String}()

    for fq in fq_

        if occursin(na, fq)

            push!(fo_, fq)

        elseif occursin(replace(na, "1" => "2"), fq)

            push!(re_, fq)

        end

    end

    n_fo, n_re = length(fo_), length(re_)

    @info "Number of forward read files = $n_fo"

    @info "Number of reverse read files = $n_re"

    if 1 >= n_fo && 1 >= n_re

        @warn "Nothing to concatenate. The number of forward reads and reverse reads are both less than or equal to 1."

    else

        dr_ = ((fo_, "R1.fastq.gz"), (re_, "R2.fastq.gz"))

        for (fq_, na) in dr_

            run(pipeline(`cat $fq_`, stdout = joinpath(pa, "$na")))

        end

        @info "Concatenated files saved at: $pa"

    end

    nothing

end

function trim(pa, r1, r2, n_jo)

    FASTQ.Support.log_sub_level_function()

    ht = joinpath(pa, "fastp.html")

    js = joinpath(pa, "fastp.json")

    ou1 = joinpath(pa, "trimmed.R1.fastq.gz")

    ou2 = joinpath(pa, "trimmed.R2.fastq.gz")

    run(
        `fastp --detect_adapter_for_pe --thread $n_jo --json $js --html $ht --in1 $r1 --in2 $r2 --out1 $ou1 --out2 $ou2`,
    )

    ou1, ou2

end

function align_dna(pa, sa, r1, r2, ge, n_jo, me)

    FASTQ.Support.log_sub_level_function()

    gei = "$ge.mmi"

    if !ispath(gei)

        run(`minimap2 -t $n_jo -d $gei $ge`)

    end

    tm = joinpath(pa, "samtools_sort")

    du = joinpath(pa, "$sa.unmarked_duplicates.bam")

    run(
        pipeline(
            `minimap2 -ax sr -t $n_jo -K $(me)G -R "@RG\tID:$sa\tSM:$sa" $gei $r1 $r2`,
            `samtools fixmate --threads $n_jo -u -m - -`,
            `samtools sort --threads $n_jo -T $tm -o $du`,
        ),
    )

    ba = joinpath(pa, "$sa.bam")

    run(`samtools markdup --threads $n_jo --reference $ge --output-fmt BAM $du $ba`)

    run(`samtools index -@ $n_jo $ba`)

    run(pipeline(`samtools stats --threads $n_jo $ba`, "$ba.stat"))

    rm(du)

    ba

end

function align_bulk_cdna_to_transcriptome(pa, r1, r2, fr, sd, tr, n_jo)

    FASTQ.Support.log_sub_level_function()

    id = "$tr.kallisto_index"

    if !ispath(id)

        @info "Creating kallisto index at $id"

        run(`kallisto index --index $id $tr`)

    end

    fu_ = ["kallisto", "quant"]

    ru_ = ["--threads", "$n_jo", "--index", "$id", "--output-dir", "$pa"]

    if r2 !== nothing

        run(`$fu_ $ru_ $r1 $r2`)

    else

        run(`$fu_ --single --fragment-length $fr --sd $sd $ru_ $r1`)

    end

end

function _index_and_stat(pa, n_jo)

    ba = joinpath(pa, "Aligned.sortedByCoord.out.bam")

    run(`samtools index -@ $n_jo $ba`)

    run(pipeline(`samtools stats --threads $n_jo $ba`, "$ba.stat"))

    ba

end

function align_bulk_cdna_to_genome(pa, r1, r2, id, n_jo)

    FASTQ.Support.log_sub_level_function()

    run(
        `star --runThreadN $n_jo --genomeDir $id --readFilesIn $r1 $r2 --readFilesCommand "gzip --decompress --stdout" --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $pa/`,
    )

    ba = _index_and_stat(pa, n_jo)

end

function align_and_quantify_bulk_cdna_to_genome(pa, r1, r2, id, n_jo)

    FASTQ.Support.log_sub_level_function()

    run(
        `star --quantMode TranscriptomeSAM GeneCounts --runThreadN $n_jo --genomeDir $id --readFilesIn $r1 $r2 --readFilesCommand "gzip --decompress --stdout" --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $pa/`,
    )

    ba = _index_and_stat(pa, n_jo)

end

function align_single_cell_cdna_to_genome(
    pas,
    r1,
    r2,
    id,
    n_jo;
    wh = nothing,
    bas = 1,
    bal = 16,
    ums = 17,
    uml = 12,
    rel = 151,
)

    FASTQ.Support.log_sub_level_function()

    if wh == nothing

        wh = joinpath(FASTQ.PR, "data", "CellRangerBarcodes", "3M-february-2018.txt")

    end

    co = "gzip --decompress --stdout"

    run(
        `star --runThreadN $n_jo --genomeDir $id --readFilesIn $r2 $r1 --readFilesCommand $co --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $pas/ --soloType Droplet --soloCBwhitelist $wh --soloCBstart $bas --soloCBlen $bal --soloUMIstart $ums --soloUMIlen $uml --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --soloBarcodeReadLength $rel`,
    )

end

end
