import gzip
import argparse
import sys

def read_vcf_metadata(vcf_file):
    metadata = {}
    opener = gzip.open if vcf_file.endswith('.gz') else open
    
    with opener(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('##'):
                key, value = line[2:].strip().split('=', 1)
                if key not in metadata:
                    metadata[key] = []
                metadata[key].append(value)
            elif line.startswith('#CHROM'):
                # This is the header line with column names
                metadata['columns'] = line[1:].strip().split('\t')
                break
            else:
                # We've reached the data lines, so we can stop
                break
    return metadata

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
                return readable_genotype
    
    return None

def print_metadata(metadata):
    print("VCF Metadata:")
    for key, values in metadata.items():
        if key == 'columns':
            print(f"  Columns: {', '.join(values)}")
        else:
            print(f"  {key}:")
            for value in values:
                print(f"    {value}")

def main():
    parser = argparse.ArgumentParser(description='Read genotype and metadata from a VCF file.')
    parser.add_argument('vcf_file', help='Path to the VCF file')
    parser.add_argument('--snp', default='rs12913832', help='Target SNP (default: rs12913832)')
    
    args = parser.parse_args()

    metadata = read_vcf_metadata(args.vcf_file)
    print_metadata(metadata)

    target_snp = args.snp
    genotype = read_genotype(args.vcf_file, target_snp)

    if genotype:
        print(f"\nYour genotype at {target_snp} is: {genotype}")
        if genotype == 'GG':
            print("This genotype is most common in people with blue eyes.")
        elif genotype in ['AA', 'GA']:
            print("This genotype is most common in people with brown eyes.")
    else:
        print(f"\nSNP {target_snp} not found in the VCF file.")

if __name__ == "__main__":
    main()