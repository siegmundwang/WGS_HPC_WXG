
STAR --runThreadN 88 --runMode genomeGenerate \
        --genomeFastaFiles <(zcat GRCh38.p13.genome.fa.gz) \
        --sjdbGTFfile <(zcat gencode.v33.annotation.gtf) \
        --sjdbOverhang 99 --genomeDir ./STAR_INDEX_GRCh38/

rsem-prepare-reference -p88 --gtf <(zcat gencode.v33.annotation.gtf) <(zcat GRCh38.p13.genome.fa.gz) ./RSEM_INDEX_GRCh38/GRCh38
