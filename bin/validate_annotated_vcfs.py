import argparse
from leiden import vcf, validation

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Outputs a VCF with all validated variants from VCFs in file_list. '
                                                 'Variants that do not validate are output to a different '
                                                 'VCF file.')

    group = parser.add_argument_group()
    group.add_argument('-f', '--file_names', required=True, help='File containing full paths to the VCF files to be processed.')
    group.add_argument('-o', '--output_file', default='lovd_validated_variants.vcf', help='Output file for validated variants (VCF).')
    group.add_argument('-d', '--discordant_output_file', default='lovd_discordant_variants.vcf', help='Output file for discordant variants (VCF).')
    args = parser.parse_args()

    total_mutation_count = 0
    total_concordant_mutation_count = 0

    vcf_header = []
    concordant_mutations = []
    discordant_mutations = []

    with open(args.file_names, 'r') as file_list:
        files_to_process = file_list.read().splitlines()

    for file in files_to_process:

        gene_mutation_count = 0
        gene_concordant_mutation_count = 0

        with open(file, 'r') as f:
            vcf_file = vcf.VCFReader(f)

            for variant in vcf_file:
                total_mutation_count += 1
                gene_mutation_count += 1

                chromosome_number = variant['CHROM']
                coordinate = variant['POS']
                viewing_interval = 25
                if coordinate != '.':
                    ucsc_link = validation.get_ucsc_location_link(chromosome_number,
                                                       str(int(coordinate) - viewing_interval),
                                                       str(int(coordinate) + viewing_interval))

                concordant_mutation_found = False
                for transcript in variant['INFO']['CSQ']:
                    if (not concordant_mutation_found) and validation.is_concordant(variant['INFO']['LOVD'][0]['PROTEIN_CHANGE'], transcript['HGVSP']):
                        total_concordant_mutation_count += 1
                        gene_concordant_mutation_count += 1
                        concordant_mutations.append(str(variant))
                        concordant_mutation_found = True

                if not concordant_mutation_found:
                    discordant_mutations.append(str(variant))

            if gene_mutation_count > 0:
                print file, gene_concordant_mutation_count, '/', gene_mutation_count, 'Concordant'
            else:
                print file, ': No annotated variants - variants could not be remapped.'

    print '-------------------------------------------'
    print total_concordant_mutation_count, '/', total_mutation_count, 'Concordant'
    print 'Concordant variants written to: ', args.discordant_output_file
    print 'Discordant variants written to: ', args.output_file

    with open(args.discordant_output_file, 'w') as discordant_file:
        discordant_file.write('\n'.join(vcf_file.header_lines) + '\n')
        discordant_file.write('\n'.join(discordant_mutations))

    with open(args.output_file, 'w') as concordant_file:
        concordant_file.write('\n'.join(vcf_file.header_lines) + '\n')
        concordant_file.write('\n'.join(concordant_mutations))