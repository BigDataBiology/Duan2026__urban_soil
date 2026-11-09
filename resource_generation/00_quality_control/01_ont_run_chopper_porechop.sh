#Use chopper and porechop to trim and remove adapters for long reads
#!/bin/bash

for n in {1..30}
  do
    gunzip -c /data/Public/soil_data/ont/sample${n}/fastq_pass.gz >sample${n}_fastq_pass.fastq

    less sample${n}_fastq_pass.fastq | chopper -q 10 -l 500 --threads 20 --contam sample${n}_ont_trim_filter_500.fastq

    porechop_abi -abi -i sample${n}_ont_trim_filter_500.fastq --format fastq --verbosity 1 --discard_middle -t 40 -o sample${n}_ont_trim_filter_500_porechop.fastq > sample${n}_ont_trim_filter_500_porechop.txt

    rm sample${n}_fastq_pass.fastq
  done