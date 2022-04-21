`number_of_jobs` _Int_<br>
Number of jobs

`memory` _Int_<br>
Gigabytes of memory

`organism` _String_<br>
human or mouse

`molecule` _String_<br>
dna or cdna

`exome` _Bool_<br>
true or false

`fragment_length` _Int_<br>
The most common cDNA fragment length in your data, used by `apply_cdna_to_genome`

`fragment_length_standard_deviation` _Int_<br>
For cdna alignment, the estimated standard deviation of fragment lengths

`sample` _String_<br>
Sample name used for file naming and name of sample column in VCF

`annotate_with_rsid` _Bool_<br>
true or false

`tool_directory` _String_<br>
Path to directory storing strelka and manta

`output_directory` _String_<br>
Path to where you want output files

`read_name_scheme` _String_<br>
Naming scheme your reads use, R1 or read1 for example

`read1` _String_<br>
Path to DNA forward read used by `apply_germline_dna_to_genome` and `apply_somatic_dna_to_genome`


`read2": "/Users/kate/craft/tool/Fastq.jl/test/data/dna/test_dna_40k.R2.fastq.gz",

`somatic_read1": "/Users/kate/craft/tool/Fastq.jl/test/data/dna/test_dna_4k.R1.fastq.gz",

`somatic_read2": "/Users/kate/craft/tool/Fastq.jl/test/data/dna/test_dna_4k.R2.fastq.gz",

`dna_read_directory": "/Users/kate/craft/tool/Fastq.jl/test/data/dna/",

`cdna_read_directory": "/Users/kate/craft/tool/Fastq.jl/test/data/cdna/",

`reference_genome": "/Users/kate/craft/data/grch/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz",

`reference_transcriptome": "/Users/kate/craft/data/mouse_reference_transcriptome/Mus_musculus.GRCm38.cdna.all.fa.gz",

`mouse_transcript_to_mouse_gene": "/Users/kate/craft/data/mouse_transcript_mouse_gene.tsv",

`chromosome_position": "/Users/kate/craft/data/grch/chromosome/chromosome.bed.gz",

`chromosome_name": "/Users/kate/craft/data/grch/chromosome/chrn_n.tsv",

`variant_database": "/Users/kate/craft/guardiome/tool/ensembl/homo_sapiens-chr1_y.vcf.gz",

`snpeff": "/Users/kate/craft/tool/Fastq.jl/tool/snpEff/snpEff.jar"
