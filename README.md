# Python-for-Tracing-Genetic-Footprints-of-the-European-Hake-
import pandas as pd
import vcfpy
import matplotlib.pyplot as plt

def extract_dna_info(vcf_file):
    # Create a VCF reader
    vcf_reader = vcfpy.Reader.from_path(vcf_file)

    # Create a list to store the DNA information
    dna_info = []

    # Iterate over each record in the VCF file
    for record in vcf_reader:
        chrom = record.CHROM
        pos = record.POS
        ref = record.REF
        alts = [alt.value for alt in record.ALT]
        info = record.INFO

        # Append the extracted information to the list
        dna_info.append({
            'Chromosome': chrom,
            'Position': pos,
            'Reference': ref,
            'Alternates': alts,
            'Info': info
        })

    # Convert the list to a DataFrame for better visualization
    df = pd.DataFrame(dna_info)
    return df

def plot_snp_distribution(df, title):
    plt.figure(figsize=(10, 6))
    df['Chromosome'].value_counts().sort_index().plot(kind='bar')
    plt.title(title)
    plt.xlabel('Chromosome')
    plt.ylabel('Number of SNPs')
    plt.show()

def main():
    vcf_file = 'example.vcf'  # Replace with your VCF file path

    # Extract DNA information
    df = extract_dna_info(vcf_file)

    # Print the extracted information
    print(df.head())

    # Plot SNP distribution across chromosomes
    plot_snp_distribution(df, 'SNP Distribution across Chromosomes')

if __name__ == "__main__":
    main()
