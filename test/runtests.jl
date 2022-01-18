using Revise

using Fastq

Fastq.test()

di = "../test/data"

re_ = Fastq.find(di)

Fastq.concatenate(re_, "R1")
