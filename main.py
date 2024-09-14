import gzip
import argparse
from Bio import SeqIO

def read_genotype(vcf_file, target_snp):
    opener = gzip.open if vcf_file.endswith('.gz') else open
    
    with opener(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            chrom, pos, id_, ref, alt, qual, filter_, info, format_ = fields[:9]
            
            if id_ == target_snp:
                genotype = fields[9].split(':')[0]
                alleles = [ref] + alt.split(',')
                readable_genotype = ''.join(alleles[int(i)] for i in genotype.replace('|', '/').split('/'))
                return chrom, int(pos), readable_genotype
    
    return None, None, None

def check_fasta_reference(fasta_file, chrom, pos):
    with gzip.open(fasta_file, 'rt') if fasta_file.endswith('.gz') else open(fasta_file, 'rt') as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id == chrom:
                return record.seq[pos-1].upper()
    return None

def main():
    parser = argparse.ArgumentParser(description='Read genotype from a VCF file for a specific SNP and check FASTA reference if not found.')
    parser.add_argument('vcf_file', help='Path to the VCF file')
    parser.add_argument('fasta_file', help='Path to the FASTA file')
    parser.add_argument('--snp', default='rs12913832', help='Target SNP (default: rs12913832)')
    parser.add_argument('--chrom', default='15', help='Chromosome number (default: 15)')
    parser.add_argument('--pos', type=int, default=28120472, help='Position in the chromosome (default: 28120472 for rs12913832 in GRCh38/hg38)')
    
    args = parser.parse_args()

    chrom, pos, genotype = read_genotype(args.vcf_file, args.snp)

    if genotype:
        print(f"Your genotype at {args.snp} is: {genotype}")
        if genotype == 'GG':
            print("This genotype is most common in people with blue eyes.")
        elif genotype in ['AA', 'GA']:
            print("This genotype is most common in people with brown eyes.")
    else:
        print(f"SNP {args.snp} not found in the VCF file. Checking FASTA reference...")
        
        # If chromosome and position weren't found in VCF, use the provided defaults
        chrom = chrom or f"chr{args.chrom}"
        pos = pos or args.pos
        
        ref_base = check_fasta_reference(args.fasta_file, chrom, pos)
        if ref_base:
            print(f"Reference base at chromosome {chrom}, position {pos}: {ref_base}")
        else:
            print(f"Unable to find the specified position in the FASTA file.")

if __name__ == "__main__":
    main()