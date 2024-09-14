import gzip
import argparse
from Bio import SeqIO

# Mapping of chromosome numbers to NCBI identifiers
CHROM_TO_NCBI = {
    '1': 'NC_000001', '2': 'NC_000002', '3': 'NC_000003', '4': 'NC_000004', '5': 'NC_000005',
    '6': 'NC_000006', '7': 'NC_000007', '8': 'NC_000008', '9': 'NC_000009', '10': 'NC_000010',
    '11': 'NC_000011', '12': 'NC_000012', '13': 'NC_000013', '14': 'NC_000014', '15': 'NC_000015',
    '16': 'NC_000016', '17': 'NC_000017', '18': 'NC_000018', '19': 'NC_000019', '20': 'NC_000020',
    '21': 'NC_000021', '22': 'NC_000022', 'X': 'NC_000023', 'Y': 'NC_000024', 'MT': 'NC_012920'
}

# Dictionary of supported SNPs with their details
SUPPORTED_SNPS = {
    'rs12913832': {
        'name': 'Eye color (HERC2 gene)',
        'chrom': '15',
        'pos': 28120472,
        'description': 'Associated with blue/brown eye color'
    },
    'rs4988235': {
        'name': 'Lactase persistence (MCM6 gene)',
        'chrom': '2',
        'pos': 136608646,
        'description': 'Associated with ability to digest lactose in adulthood'
    }
}

def read_genotype_vcf(vcf_file, target_snp):
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
    ncbi_chrom = CHROM_TO_NCBI.get(chrom.replace('chr', ''))
    if not ncbi_chrom:
        print(f"Unable to map chromosome {chrom} to NCBI identifier.")
        return None

    with gzip.open(fasta_file, 'rt') if fasta_file.endswith('.gz') else open(fasta_file, 'rt') as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id.startswith(ncbi_chrom):
                return record.seq[pos-1].upper()
    return None

def interpret_result(snp, vcf_genotype, ref_base):
    if snp == 'rs12913832':
        if vcf_genotype:
            print(f"Your genotype at {snp} is: {vcf_genotype}")
            if vcf_genotype == 'GG':
                print("YOU LIKELY HAVE BLUE EYES.")
                print("This GG genotype is strongly associated with blue eye color.")
            elif vcf_genotype in ['AA', 'GA']:
                print("YOU LIKELY HAVE BROWN EYES.")
                print("This AA or GA genotype is strongly associated with brown eye color.")
        elif ref_base:
            print(f"SNP {snp} was not found in your VCF file. The reference base is '{ref_base}'.")
            if ref_base == 'A':
                print("YOUR GENOTYPE IS MOST LIKELY AA.")
                print("This means you probably have brown eyes.")
                print("However, if your VCF has poor coverage in this region, your actual genotype could be different.")
            elif ref_base == 'G':
                print("YOUR GENOTYPE IS POSSIBLY GG, BUT THIS IS LESS CERTAIN.")
                print("If GG, you might have blue eyes, but this is less common for the reference genome.")
                print("Consider getting this region sequenced with better coverage for a more definitive answer.")
        else:
            print(f"Unable to find information for {snp} in both VCF and reference genome.")
    
    elif snp == 'rs4988235':
        if vcf_genotype:
            print(f"Your genotype at {snp} is: {vcf_genotype}")
            if vcf_genotype in ['AA', 'GA']:
                print("YOU LIKELY HAVE LACTASE PERSISTENCE.")
                print("This means you can probably digest lactose (milk sugar) as an adult.")
                if vcf_genotype == 'GA':
                    print("Even with one copy of the 'A' allele, you're likely to have lactase persistence.")
                else:  # AA
                    print("With two copies of the 'A' allele, you have a very high likelihood of lactase persistence.")
            elif vcf_genotype == 'GG':
                print("YOU ARE LIKELY LACTOSE INTOLERANT.")
                print("This means you may have difficulty digesting lactose (milk sugar) as an adult.")
        elif ref_base:
            print(f"SNP {snp} was not found in your VCF file. The reference base is '{ref_base}'.")
            if ref_base == 'G':
                print("YOUR GENOTYPE IS MOST LIKELY GG.")
                print("This means you are PROBABLY LACTOSE INTOLERANT.")
                print("However, if your VCF has poor coverage in this region, your actual genotype could be different.")
                print("If you can digest dairy products without issues, you might actually have a GA or AA genotype.")
            elif ref_base == 'A':
                print("YOUR GENOTYPE IS POSSIBLY AA OR GA, BUT THIS IS LESS CERTAIN.")
                print("If AA or GA, you likely have lactase persistence, but this is less common for the reference genome.")
                print("Consider getting this region sequenced with better coverage for a more definitive answer.")
        else:
            print(f"Unable to find information for {snp} in both VCF and reference genome.")
    
    
    else:
        if vcf_genotype:
            print(f"Your genotype at {snp} is: {vcf_genotype}")
            print("No specific interpretation is available for this SNP.")
        elif ref_base:
            print(f"SNP {snp} was not found in your VCF file. The reference base is '{ref_base}'.")
            print("No specific interpretation is available for this SNP.")
        else:
            print(f"Unable to find information for {snp} in both VCF and reference genome.")

def list_supported_snps():
    print("Supported SNPs:")
    for i, (snp, info) in enumerate(SUPPORTED_SNPS.items(), 1):
        print(f"{i}. {snp} - {info['name']}: {info['description']}")

def get_user_snp_choice():
    while True:
        choice = input("Enter the number of the SNP you want to analyze (or 'q' to quit): ")
        if choice.lower() == 'q':
            return None
        try:
            choice = int(choice)
            if 1 <= choice <= len(SUPPORTED_SNPS):
                return list(SUPPORTED_SNPS.keys())[choice - 1]
        except ValueError:
            pass
        print("Invalid choice. Please try again.")

def main():
    parser = argparse.ArgumentParser(description='Analyze genotype from VCF and FASTA files')
    parser.add_argument('--vcf', help='Path to the VCF file')
    parser.add_argument('--fasta', help='Path to the FASTA file')
    parser.add_argument('--snp', help='Target SNP (e.g., rs12913832)')
    parser.add_argument('--chrom', help='Chromosome number')
    parser.add_argument('--pos', type=int, help='Position in the chromosome')
    
    args = parser.parse_args()

    if not args.snp:
        list_supported_snps()
        chosen_snp = get_user_snp_choice()
        if not chosen_snp:
            print("Exiting program.")
            return
        args.snp = chosen_snp
        args.chrom = SUPPORTED_SNPS[chosen_snp]['chrom']
        args.pos = SUPPORTED_SNPS[chosen_snp]['pos']

    vcf_genotype = None
    ref_base = None

    if args.vcf:
        chrom, pos, vcf_genotype = read_genotype_vcf(args.vcf, args.snp)
        if not vcf_genotype:
            print(f"SNP {args.snp} not found in the VCF file.")
        else:
            print(f"SNP {args.snp} not found in the VCF file! Value: {vcf_genotype}")
    
    if args.fasta:
        ref_base = check_fasta_reference(args.fasta, args.chrom, args.pos)
        if ref_base:
            print(f"\nFASTA Analysis: Reference base at chromosome {args.chrom}, position {args.pos}: {ref_base}")
        else:
            print(f"\nUnable to find the specified position in the FASTA file.")
    
    interpret_result(args.snp, vcf_genotype, ref_base)

    if not (args.vcf or args.fasta):
        print("Please provide at least one input file (VCF or FASTA).")

if __name__ == "__main__":
    main()