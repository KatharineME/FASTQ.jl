`number_of_jobs`<br>
Int, number of jobs

#### `memory`

Int, Gigabytes of memory

#####`organism`
String, human or mouse

"molecule": String, dna or cdna
"exome": Bool, true or false
"fragment_length": Int, For cdna alignment to genome, the most common fragment length
"fragment_length_standard_deviation": Int, For cdna alignment, the estimated standard deviation of fragment lengths
"sample": (String) sample name used for file naming and name of sample column in VCF
"annotate_with_rsid": (Bool) true or false
"tool_directory": (String) Path to directory storing strelka and manta
"output_directory": (String) Path to where you want output files
"read_name_scheme": (String) Naming scheme your reads use, R1 or read1 for example
"read1": "/Users/kate/craft/tool/Fastq.jl/test/data/dna/test_dna_40k.R1.fastq.gz",
"read2": "/Users/kate/craft/tool/Fastq.jl/test/data/dna/test_dna_40k.R2.fastq.gz",
"somatic_read1": "/Users/kate/craft/tool/Fastq.jl/test/data/dna/test_dna_4k.R1.fastq.gz",
"somatic_read2": "/Users/kate/craft/tool/Fastq.jl/test/data/dna/test_dna_4k.R2.fastq.gz",
"dna_read_directory": "/Users/kate/craft/tool/Fastq.jl/test/data/dna/",
"cdna_read_directory": "/Users/kate/craft/tool/Fastq.jl/test/data/cdna/",
"reference_genome": "/Users/kate/craft/data/grch/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz",
"reference_transcriptome": "/Users/kate/craft/data/mouse_reference_transcriptome/Mus_musculus.GRCm38.cdna.all.fa.gz",
"mouse_transcript_to_mouse_gene": "/Users/kate/craft/data/mouse_transcript_mouse_gene.tsv",
"chromosome_position": "/Users/kate/craft/data/grch/chromosome/chromosome.bed.gz",
"chromosome_name": "/Users/kate/craft/data/grch/chromosome/chrn_n.tsv",
"variant_database": "/Users/kate/craft/guardiome/tool/ensembl/homo_sapiens-chr1_y.vcf.gz",
"snpeff": "/Users/kate/craft/tool/Fastq.jl/tool/snpEff/snpEff.jar"
