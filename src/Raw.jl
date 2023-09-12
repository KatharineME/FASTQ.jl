module Raw

using FASTQ

function find(di)

    FASTQ.Support.log_sub_level_function()

    fq_ = Vector{String}()

    na_n = Dict(".fq" => 0, ".fastq" => 0, "fq.gz" => 0, "fastq.gz" => 0)

    for (ro, di_, fi_) in walkdir(di)

        for fi in fi_

            for (na, _) in na_n

                if endswith(fi, na)

                    na_n[na] += 1

                    if endswith(fi, ".gz")

                        push!(fq_, joinpath(ro, fi))

                    end

                end

            end

        end

    end

    @info "File types found in $di"

    for na in na_n

        @info na

    end

    @info "Size of gzipped files"

    for fi in fq_

        @info "File $fi is:" Base.format_bytes(stat(fi).size)

    end

    fq_

end

function check_read(pa, fq_, n_jo)

    FASTQ.Support.log_sub_level_function()

    FASTQ.Support.trash_remake_directory(pa)

    th = minimum((length(fq_), n_jo))

    run(`fastqc --threads $th --outdir $pa $fq_`)

    run(`multiqc --outdir $pa $pa`)

end

function check_read(pa, r1, r2, n_jo; sor1 = nothing, sor2 = nothing)

    if sor1 === nothing

        fq_ = [r1, r2]

    else

        fq_ = [r1, r2, sor1, sor2]

    end

    check_read(joinpath(pa, "check_raw"), fq_, n_jo)

end

function concatenate(fq_; na = "R1")

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

    n_fo = length(fo_)

    n_re = length(re_)

    @info "Number of forward read files = $n_fo"

    @info "Number of reverse read files = $n_re"

    fo1 = split(fo_[1], "/")

    co = ""

    for i in reverse(1:length(fo1))
        if all(split(fq, "/")[i] == fo1[i] for fq in fq_)
            co = fo1[i]
            break
        end
    end

    pac = joinpath(split(fo_[1], co)[1], "$(co)Concatenated")

    if 1 >= n_fo && 1 >= n_re

        @warn "Nothing to concatenate. The umber of forward reads and reverse reads are both less than or equal to 1."

    else

        FASTQ.Support.trash_remake_directory(pac)

        @info "Concatenating"

        dr_ = ((fo_, "R1.fastq.gz"), (re_, "R2.fastq.gz"))

        for (fq_, na) in dr_

            run(pipeline(`cat $fq_`, stdout = joinpath(pac, "$na")))

        end

        @info "Concatenated files saved at: $pac"

    end

end

function trim(pa, n_jo, r1, r2)

    FASTQ.Support.log_sub_level_function()

    FASTQ.Support.trash_remake_directory(pa)

    ht = joinpath(pa, "fastp.html")

    js = replace(ht, "html" => "json")

    ou1 = joinpath(pa, FASTQ._TR1)

    ou2 = joinpath(pa, FASTQ._TR2)

    run(
        `fastp --detect_adapter_for_pe --thread $n_jo --json $js --html $ht --in1 $r1 --in2 $r2 --out1 $ou1 --out2 $ou2`,
    )

end

function psuedoalign(ou, tr, n_jo, r1, r2, fr, sd)

    FASTQ.Support.log_sub_level_function()

    id = "$tr.kallisto_index"

    if !ispath(id)

        @info "Creating kallisto index at $id"

        run(`kallisto index --index $id $tr`)

    end

    fu_ = ["kallisto", "quant"]

    ru_ = ["--threads", "$n_jo", "--index", "$id", "--output-dir", "$ou"]

    if r2 !== nothing

        run(`$fu_ $ru_ $r1 $r2`)

    else

        run(`$fu_ --single --fragment-length $fr --sd $sd $ru_ $r1`)

    end

end

function align_cdna(pa, ge, n_jo, sa, r1, r2)

    FASTQ.Support.log_sub_level_function()

    FASTQ.Support.trash_remake_directory(pa)

    id = joinpath(dirname(ge), "star_indexes")

    if !ispath(id)

        mkdir(id)

        @info "Making STAR indices"

        ged = splitext(ge)[1]

        if !isfile(ged)

            run(pipeline(`bgzip --decompress --stdout $ge`, `$ged`))

        end

        run(
            `star --runThreadN $n_jo --runMode genomeGenerate --genomeDir $id --genomeFastaFiles $ged`,
        )

    end

    pr = joinpath(pa, "$(sa).")

    run(
        `star --runThreadN $n_jo --genomeDir $id --readFilesIn $r1 $r2 --readFilesCommand "gzip --decompress --stdout" --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $pr`,
    )

    ba = "$(pr)Aligned.sortedByCoord.out.bam"

    run(`samtools index -@ $n_jo $ba`)

    run(pipeline(`samtools stats --threads $n_jo $ba`, "$ba.stat"))

end

function align_cdna(pa, cd, re, n_jo; al = "transcriptome", fr = 51, sd = 0.05)

    FASTQ.Support.log_sub_level_function()

    fq_ = find(cd)

    na_ = ["R1", "read1", "_1.fq"]

    for fq1 in fq_

        for na in na_

            if occursin(na, fq1)

                fq2 = replace(fq1, na => replace(na, "1" => "2"))

                if !isfile(fq2)

                    fq2 = nothing

                end

                sa = last(splitdir(splitext(split(fq1, na)[1])[1]))

                @info "Working on sample: $sa"

                pas = joinpath(pa, sa)

                if al == "transcriptome"

                    psuedoalign(pas, re, n_jo, fq1, fq2, fr, sd)

                elseif al == "genome"

                    align_cdna(pas, re, n_jo, sa, fq1, fq2)

                end

            end

        end

    end

end


function align_dna(pa, sa, ba, r1, r2, ge, n_jo, me)

    FASTQ.Support.log_sub_level_function()

    FASTQ.Support.trash_remake_directory(pa)

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

    run(`samtools markdup --threads $n_jo --reference $ge --output-fmt BAM $du $ba`)

    run(`samtools index -@ $n_jo $ba`)

    run(pipeline(`samtools stats --threads $n_jo $ba`, "$ba.stat"))

    rm(du)

end

function align_single_cell_cdna(pa, sa, r1, r2, ge, n_jo)

    FASTQ.Support.log_sub_level_function()

    FASTQ.Support.trash_remake_directory(pa)

    # Change to star index made with gtf gene annotation file
    id = joinpath(dirname(ge), "star_indexes")

    if !ispath(id)

        mkdir(id)

        ged = splitext(ge)[1]

        if !isfile(ged)

            run(pipeline(`bgzip --decompress --stdout $ge`, `$ged`))

        end

        run(
            `star --runThreadN $n_jo --runMode genomeGenerate --genomeDir $id --genomeFastaFiles $ged`,
        )

    end

    pr = joinpath(pa, "$(sa).")

    # Update star run
    run(
        `star --runThreadN $n_jo --genomeDir $id --readFilesIn $r1 $r2 --readFilesCommand "gzip --decompress --stdout" --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $pr`,
    )

    ba = "$(pr)Aligned.sortedByCoord.out.bam"

    run(`samtools index -@ $n_jo $ba`)

    run(pipeline(`samtools stats --threads $n_jo $ba`, "$ba.stat"))

end

end
