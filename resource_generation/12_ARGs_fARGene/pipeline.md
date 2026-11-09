### Split mag sequences in several files
```bash
cd /work/microbiome/users/juan/UrbanSoil/Genes/mags/data
python ../split_data.py /work/microbiome/UrbanSoil/data/UrbanSoilGenes/gene_prediction/mags/mags.faa.gz 10000000 &
```

### Run fARGene in contigs data 
#### Run fARGene
```bash
cd /work/microbiome/users/juan/UrbanSoil/contigs/
conda activate snakemake-env
snakemake -s Snakefile -p -j 30 --use-conda --profile aqua
```

#### Compile the results
```bash
cd /work/microbiome/users/juan/UrbanSoil/contigs/
for f in fargene/*@@@*/predictedGenes/predicted-orfs-amino.fasta; do
    sample_class=$(echo "$f" | cut -d/ -f2)
    sample=${sample_class%%@@@*}
    class=${sample_class##*@@@}

    sample_short=${sample%%_*}
    prefix="${sample_short}@@@${class}"

    awk -v prefix="$prefix" '
    /^>/{
        gsub(/^>/,"")
        gsub(/orfs-translated_/,"")

        n=split($0,a,"_")
        new=a[1]
        for(i=2;i<=n-2;i++){
            new=new"_"a[i]
        }

        print ">"prefix"@@@"new
        next
    }
    {print}
    ' "$f" >> predicted-orfs-amino.fasta
done


for f in fargene/*@@@*/predictedGenes/predicted-orfs.fasta; do
    sample_class=$(echo "$f" | cut -d/ -f2)
    sample=${sample_class%%@@@*}
    class=${sample_class##*@@@}

    sample_short=${sample%%_*}
    prefix="${sample_short}@@@${class}"

        awk -v prefix="$prefix" '
    /^>/{
        gsub(/^>/,"")

        sub(/^.*long-orfs_/,"")
        sub(/_seq.*/,"")    

        print ">"prefix"@@@"$0
        next
    }
    {print}
    
    ' "$f" >> predicted-orfs.fasta

done

for f in fargene/*@@@*/predictedGenes/*-filtered-peptides.fasta; do

    sample_class=$(echo "$f" | cut -d/ -f2)
    sample=${sample_class%%@@@*}
    class1=${sample_class##*@@@}

    sample_short=${sample%%_*}
    prefix="${sample_short}@@@${class1}"

    awk -v prefix="$prefix" '
    /^>/{
        gsub(/^>/,"")

        sub(/^.*\.gz/,"")   # remove everything up to .gz
        sub(/_seq.*/,"")    # remove _seq and everything after

        print ">"prefix"@@@"$0
        next
    }
    {print}
    ' "$f" >> filter-contigs-amino.fasta

done

for f in fargene/*@@@*/predictedGenes/*-filtered.fasta; do

    sample_class=$(echo "$f" | cut -d/ -f2)
    sample=${sample_class%%@@@*}
    class1=${sample_class##*@@@}

    sample_short=${sample%%_*}
    prefix="${sample_short}@@@${class1}"

    awk -v prefix="$prefix" '
    /^>/{
        gsub(/^>/,"")

        sub(/^.*\.gz/,"")   # remove everything up to .gz
        sub(/_seq.*/,"")    # remove _seq and everything after

        print ">"prefix"@@@"$0
        next
    }
    {print}
    ' "$f" >> filter-contigs.fasta

done

```



### Rename the sequences to avoid duplicate names
```bash
conda activate seqkit
seqkit rename -n predicted-orfs.fasta -o predicted-orfs-nonredundant.fasta
seqkit rename -n predicted-orfs-amino.fasta -o predicted-orfs-amino-nonredundant.fasta
seqkit rename -n filter-contigs.fasta -o filter-contigs-nonredundant.fasta
seqkit rename -n filter-contigs-amino.fasta -o filter-contigs-amino-nonredundant.fasta

```


### Run CD-HIT on the amino acid sequences, VSEARCH for nucleotide sequences
```bash
conda activate cdhit_env
cd clustering_contigs_amino
cd-hit -i ../filter-contigs-amino-nonredundant.fasta -o cluster95-contigs-amino-nonredundant -c 0.95 -n 5 -M 0 -T 8

cd ../clustering_orfs_amino
cd-hit -i ../predicted-orfs-amino-nonredundant.fasta -o cluster95-orfs-amino-nonredundant -c 0.95 -n 5 -M 0 -T 8

conda activate vsearch
cd ../clustering_orfs_amino
vsearch --cluster_fast ../predicted-orfs-nonredundant.fasta --id 0.95 --centroids cluster95-orfs-nonredundant.fasta --uc cluster95-orfs-nonredundant.uc 

cd ../clustering_contigs
vsearch --cluster_fast ../filter-contigs-nonredundant.fasta --id 0.95 --centroids cluster95-contigs-nonredundant.fasta --uc cluster95-contigs-nonredundant.uc 

```


### BLAST the predicted genes to the centroids to finish clustering with CD-HIT
```bash
cd /work/microbiome/users/juan/UrbanSoil/contigs/clustering_orfs_amino
conda activate diamond
diamond makedb --in cluster95-orfs-amino-nonredundant -d ref_db
diamond blastp -q ../predicted-orfs-amino-nonredundant.fasta -d ref_db.dmnd -o closest_hits_orfs_amino.tsv --max-target-seqs 1 --outfmt 6 qseqid sseqid pident evalue bitscore

cd /work/microbiome/users/juan/UrbanSoil/contigs/clustering_contigs_amino
conda activate diamond
diamond makedb --in cluster95-contigs-amino-nonredundant -d ref_db
diamond blastp -q ../filter-contigs-amino-nonredundant.fasta -d ref_db.dmnd -o closest_hits_orfs_amino.tsv --max-target-seqs 1 --outfmt 6 qseqid sseqid pident evalue bitscore

```

### Running ResFinder on fARGene results 
```bash

cd ../resfinder_orfs/
conda activate resfinder
export CGE_RESFINDER_RESGENE_DB="/work/microbiome/users/juan/resfinder_databases/resfinder_db"
export CGE_RESFINDER_RESPOINT_DB="/work/microbiome/users/juan/resfinder_databases/pointfinder_db/"
export CGE_DISINFINDER_DB="/work/microbiome/users/juan/resfinder_databases/disinfinder_db"

python -m resfinder -o results -l 0.6 -t 0.8 --acquired -ifa ../predicted-orfs-nonredundant.fasta
awk -F"\t" '{print $1 "\t" $2 "\t" $3  "\t" $4 "\t" $6 }' results/ResFinder_results_tab.txt > resfinfer_on_orfs-nonredundant.tsv

# same for contigs and contigs amino, and resfinder_orf_amino
```


### Running eggNOG-mapper on fARGene results
```bash
cd eggnog_orfs_amino
conda activate eggnog_out
export EGGNOG_DATA_DIR="/work/microbiome/users/juan/eggnog_data"
export PATH=/home/user/eggnog-mapper:/work/microbiome/users/juan/eggnog-mapper/eggnogmapper/bin:"$PATH"
../../../eggnog-mapper/emapper.py -i ../clustering_orfs_amino/cluster95-orfs-amino-nonredundant -o eggnog_out
```











### Run RGI on the genes from MAGs
```bash
cd /work/microbiome/users/juan/UrbanSoil/Genes/mags/rgi
zcat /work/microbiome/UrbanSoil/data/UrbanSoilGenes/gene_prediction/mags/mags.faa.gz | sed 's/\*//g' | gzip > mags_clean.faa.gz
conda activate RGI
rgi main -a DIAMOND -i mags_clean.faa.gz -o rgi  --local --clean -t protein -n 16
```


### Running AMRFinderPlus on fARGene results
```bash
cd /work/microbiome/users/juan/UrbanSoil/Genes/contigs/clusters
conda activate seqkit
seqkit translate --frame 6 centroids.fasta |seqkit rename -n |seqkit sort -l -r | seqkit replace -p "_[\d]+$" -r "" | seqkit rmdup -n > centroids_longest_orf.faa

cd ../amrfinder
conda activate amrfinder
amrfinder -p ../clusters/centroids_longest_orf.faa > amrfinder_on_fargene_centroids95.tsv

conda activate fargene
prodigal -i ../clusters/centroids.fasta -a prodigal_centroids.faa  -p meta
amrfinder -p prodigal_centroids.faa > amrfinder_on_fargene_centroids95_prodigal.tsv


cd /work/microbiome/users/juan/UrbanSoil/Genes/mags/amrfinder
amrfinder -p ../cdhit/cluster95 > amrfinder_on_fargene_centroids95.tsv

```



### Run fARGene in mags data 
```bash
cd /work/microbiome/users/juan/UrbanSoil/contigs/
snakemake -s Snakefile_mags -p -j 20 --use-conda



cd /work/microbiome/users/juan/UrbanSoil/Genes/mags
for filename in fargene_predicted_suffix/*.fasta; do
    tag="${filename##*/}"
    #tag="${tag#*.gz-}"
    #tag="${tag%.fasta}"
    awk -v tag="$tag" '
        /^>/ { print ">" substr($0, 2) "@@@" tag; next }
        { print }
    ' "$filename" >> predicted_genes_with_class_clustering.fasta
done
```

### Run CD-HIT on the predicted genes from MAGs
```bash
conda activate cdhit-env
cd /work/microbiome/users/juan/UrbanSoil/Genes/mags
cd-hit -i ../predicted_genes_with_class_clustering.fasta -o cluster95 -c 0.95 -n 5 -M 0 -T 8
```
