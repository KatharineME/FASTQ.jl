module Raw

using ..FASTQ

function find(pa)

    fi_ = readdir(pa; join = true)

    fi_[findall(occursin.(".gz", fi_))]

end

function test(pa, n_jo)

    fi_ = find(pa)

    r1 = RN1

    r2 = replace(r1, "1" => "2")

    li_ = ()

    for fi in fi_

        @info "Checking $fi"

        if occursin(r1, fi) | occursin(r2, fi)

            @info "File named correctly: $fi."

        else

            @warn "File named incorrectly, does not contain $r1 or $r2: $fi."

        end

        @info "Testing bgzip integrity"

        try
            run(`bgzip -t --threads $n_jo $fi`)

            @info "bzip test passed."

        catch er

            #wc -l

        end

    end

end

function check(pa, fq_, n_jo)

    FASTQ.Support.log_sub_level_function()

    run(`fastqc --threads $(minimum((lastindex(fq_), n_jo))) --outdir $pa $fq_`)

    run(`multiqc --outdir $pa $pa`)

    fq_[1], fq_[2]

end

function concatenate(pa, fq_; na = FASTQ._RN1)

    FASTQ.Support.log_sub_level_function()

    fo_ = String[]

    re_ = String[]

    for fq in fq_

        if occursin(na, fq)

            push!(fo_, fq)

        elseif occursin(replace(na, "1" => "2"), fq)

            push!(re_, fq)

        end

    end

    n_fo, n_re = lastindex(fo_), lastindex(re_)

    @info "Number of forward read files = $n_fo"

    @info "Number of reverse read files = $n_re"

    if 1 >= n_fo && 1 >= n_re

        @warn "Nothing to concatenate. The number of forward reads and reverse reads are both less than or equal to 1."

    else

        dr_ = ((fo_, "$(FASTQ._RN1).fastq.gz"), (re_, "$(FASTQ._RN2).fastq.gz"))

        for (fq_, na) in dr_

            run(pipeline(`cat $fq_`; stdout = joinpath(pa, "$na")))

        end

        @info "Concatenated files saved at: $pa"

    end

    nothing

end

function trim(pa, r1, r2, n_jo)

    FASTQ.Support.log_sub_level_function()

    ou1 = joinpath(pa, "trimmed.$(FASTQ._RN1).fastq.gz")

    ou2 = replace(ou1, FASTQ._RN1 => FASTQ._RN2)

    run(
        `fastp --in1 $r1 --in2 $r2 --out1 $ou1 --out2 $ou2 --detect_adapter_for_pe --cut_tail --trim_poly_g --trim_poly_x --length_required 48 --thread $n_jo --json $(joinpath(pa, "fastp.json")) --html $(joinpath(pa, "fastp.html"))`,
    )

    ou1, ou2

end

function align_dna(pa, sa, r1, r2, ge, n_jo, me)

    FASTQ.Support.log_sub_level_function()

    !ispath("$ge.bwt") ? run(`bwa index $ge`) : nothing

    tm = joinpath(pa, "samtools_sort")

    du = joinpath(pa, "$sa.unmarked_duplicates.bam")

    #TODO: Check bwa mem call
    run(
        pipeline(
            `bwa mem -t $n_jo -v 3 -R "@RG\tID:$sa\tSM:$sa" $ge $r1 $r2`,
            `samtools fixmate --threads $n_jo -u -m - -`,
            `samtools sort --threads $n_jo -T $tm -o $du`,
        ),
    )

    ba = joinpath(pa, "$sa.bam")

    run(
        `samtools markdup --threads $n_jo --reference $ge --output-fmt BAM $du $ba`,
    )

    run(`samtools index -@ $n_jo $ba`)

    run(pipeline(`samtools stats --threads $n_jo $ba`, "$ba.stat"))

    rm(du)

    ba

end

# function align_dna(pa, sa, r1, r2, ge, n_jo, me)
# 
#     FASTQ.Support.log_sub_level_function()
# 
#     gei = "$ge.mmi"
# 
#     !ispath(gei) ? run(`minimap2 -t $n_jo -d $gei $ge`) : nothing
# 
#     tm = joinpath(pa, "samtools_sort")
# 
#     du = joinpath(pa, "$sa.unmarked_duplicates.bam")
# 
#     run(
#         pipeline(
#             `minimap2 -ax sr -t $n_jo -K $(me)G -R "@RG\tID:$sa\tSM:$sa" $gei $r1 $r2`,
#             `samtools fixmate --threads $n_jo -u -m - -`,
#             `samtools sort --threads $n_jo -T $tm -o $du`,
#         ),
#     )
# 
#     ba = joinpath(pa, "$sa.bam")
# 
#     run(`samtools markdup --threads $n_jo --reference $ge --output-fmt BAM $du $ba`)
# 
#     run(`samtools index -@ $n_jo $ba`)
# 
#     run(pipeline(`samtools stats --threads $n_jo $ba`, "$ba.stat"))
# 
#     rm(du)
# 
#     ba
# 
# end

function align_bulk_cdna_to_transcriptome(pa, r1, r2, fr, sd, tr, n_jo)

    FASTQ.Support.log_sub_level_function()

    id = "$tr.kallisto_index"

    if !ispath(id)

        @info "Creating kallisto index at $id"

        run(`kallisto index --index $id $tr`)

    end

    fu_ = ("kallisto", "quant")

    ru_ = ("--threads", "$n_jo", "--index", "$id", "--output-dir", "$pa")

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

    _index_and_stat(pa, n_jo)

end

function align_and_quantify_bulk_cdna_to_genome(pa, r1, r2, id, n_jo)

    FASTQ.Support.log_sub_level_function()

    run(
        `star --quantMode TranscriptomeSAM GeneCounts --runThreadN $n_jo --genomeDir $id --readFilesIn $r1 $r2 --readFilesCommand "gzip --decompress --stdout" --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $pa/`,
    )

    _index_and_stat(pa, n_jo)

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

    if wh === nothing
        wh = joinpath(
            FASTQ.PR,
            "data",
            "CellRangerBarcodes",
            "3M-february-2018.txt",
        )
    else
        nothing
    end

    run(
        `star --runThreadN $n_jo --genomeDir $id --readFilesIn $r2 $r1 --readFilesCommand $("gzip --decompress --stdout") --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $pas/ --soloType Droplet --soloCBwhitelist $wh --soloCBstart $bas --soloCBlen $bal --soloUMIstart $ums --soloUMIlen $uml --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --soloBarcodeReadLength $rel`,
    )

end

end
