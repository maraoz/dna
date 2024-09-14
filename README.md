# DNA Tools

### Install
```sh
pip install -r requirements.txt
```

### Required data
1) Your VCF file: you should get this from a full genome sequencing company.

2) Your reference FASTA file: you can check your VCF metadata for which reference genome it uses.
In most cases it will be hg38, which you can download here:
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.fna.gz

### Run:
```sh
python3 main.py --vcf path_to_your.vcf --fasta GCF_000001405.26_GRCh38_genomic.fna
```