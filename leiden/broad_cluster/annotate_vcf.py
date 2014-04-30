import os
import subprocess

# VCF-Annotate variant data sources
MAC26K_AC = '/humgen/atgu1/fs03/lek/vcf-annotate/mac26k_AC.txt.gz'
EUR_AC = '/humgen/atgu1/fs03/lek/vcf-annotate/eur_AC.txt.gz'
SAS_AC = '/humgen/atgu1/fs03/lek/vcf-annotate/sas_AC.txt.gz'
AFR_AC = '/humgen/atgu1/fs03/lek/vcf-annotate/afr_AC.txt.gz'
AMR_AC = '/humgen/atgu1/fs03/lek/vcf-annotate/amr_AC.txt.gz'
EAS_AC = '/humgen/atgu1/fs03/lek/vcf-annotate/eas_AC.txt.gz'
HGMD = '/humgen/atgu1/fs03/lek/HGMD/vcf-annotate/hgmd2013.txt.gz'
HGMD_GENES = '/humgen/atgu1/fs03/lek/HGMD/vcf-annotate/hgmd_genes_anno.txt.gz'
DBSNP = '/humgen/atgu1/fs03/lek/vcf-annotate/dbsnp_137.b37.txt.gz'

# Information about the available VCF-Annotate data sources
MAC26K_COLUMNS = "CHROM,FROM,REF,ALT,INFO/AC_MAC26K,INFO/AN_MAC26K"
MAC26K_DESCRIPTION = "key=INFO,ID=AC_MAC26K,Number=1,Type=String,Description='26KJoint Called AC'"
MAC26K_DESCRIPTION2 = "key=INFO,ID=AN_MAC26K,Number=1,Type=String,Description='26KJoint Called AN'"

EUR_COLUMNS = "CHROM,FROM,REF,ALT,INFO/AC_EUR,INFO/AN_EUR"
EUR_DESCRIPTION = "key=INFO,ID=AC_EUR,Number=1,Type=String,Description='EuropeanAC'"
EUR_DESCRIPTION2 = "key=INFO,ID=AN_EUR,Number=1,Type=String,Description='EuropeanAN'"

SAS_COLUMNS = "CHROM,FROM,REF,ALT,INFO/AC_SAS,INFO/AN_SAS"
SAS_DESCRIPTION = "key=INFO,ID=AC_SAS,Number=1,Type=String,Description='South Asian AC'"
SAS_DESCRIPTION2 = "key=INFO,ID=AN_SAS,Number=1,Type=String,Description='South Asian AN'"

AFR_COLUMNS = "CHROM,FROM,REF,ALT,INFO/AC_AFR,INFO/AN_AFR"
AFR_DESCRIPTION = "key=INFO,ID=AC_AFR,Number=1,Type=String,Description='African AC'"
AFR_DESCRIPTION2 = "key=INFO,ID=AN_AFR,Number=1,Type=String,Description='African AN'"

AMR_COLUMNS = "CHROM,FROM,REF,ALT,INFO/AC_AMR,INFO/AN_AMR"
AMR_DESCRIPTION = "key=INFO,ID=AC_AMR,Number=1,Type=String,Description='AmericanAC'"
AMR_DESCRIPTION2 = "key=INFO,ID=AN_AMR,Number=1,Type=String,Description='AmericanAN'"

EAS_COLUMNS = "CHROM,FROM,REF,ALT,INFO/AC_EAS,INFO/AN_EAS"
EAS_DESCRIPTION = "key=INFO,ID=AC_EAS,Number=1,Type=String,Description='East Asian AC'"
EAS_DESCRIPTION2 = "key=INFO,ID=AN_EAS,Number=1,Type=String,Description='East Asian AN'"

HGMD_SITES_COLUMNS = "CHROM,FROM,-,-,INFO/HGMD_SITE"
HGMD_SITES_DESCRIPTION = "key=INFO,ID=HGMD_SITE,Number=1,Type=String,Description='HGMD SNP Site'"

HGMD_MUTATIONS_COLUMNS = "CHROM,FROM,REF,ALT,INFO/HGMD_MUT"
HGMD_MUTATIONS_DESCRIPTION = "key=INFO,ID=HGMD_MUT,Number=1,Type=String,Description='HGMD SNP Mutation'"

HGMD_GENES_COLUMNS = "CHROM,FROM,TO,INFO/HGMD_GENE,-,-,-"
HGMD_GENES_DESCRIPTION = "key=INFO,ID=HGMD_GENE,Number=0,Type=Flag,Description='HGMD Gene'"

DBSNP_COLUMNS = "CHROM,FROM,REF,ALT,INFO/DBSNP"
DBSNP_DESCRIPTION = "key=INFO,ID=DBSNP,Number=1,Type=String,Description='DBSNP'"


def annotate_vep(input_file, output_file):
    """
    Annotate VCF file with Variant Effect Predictor.

    Args:
        input_file (str): input VCF file path
        output_file (str): output VCF file path (VEP annotation added to file).

    """

    vep_path = '/humgen/atgu1/fs03/lek/VEP25/'
    vep_exe_path = os.path.join(vep_path, 'variant_effect_predictor.pl')

    subprocess.call(['perl',
                     vep_exe_path,
                     '--offline',
                     '--dir', vep_path,
                     '--format', 'vcf',
                     '--vcf',
                     '--geuvadis',
                     '--force_overwrite',
                     '--no_progress',
                     '--input_file', input_file,
                     '--output_file', output_file])


def annotate_26k(input_file, output_file):
    """
    Annotate VCF file with MacArthur 26K data.

    Args:
        input_file (str): input VCF file path
        output_file (str): output VCF file path (26K allele frequency annotation added to file)

    """

    input_text = subprocess.Popen(['cat', input_file], stdout=subprocess.PIPE)

    annotate_mac26k = vcf_annotate(input_text, MAC26K_AC, MAC26K_COLUMNS, [MAC26K_DESCRIPTION, MAC26K_DESCRIPTION2])
    annotate_eur = vcf_annotate(annotate_mac26k, EUR_AC, EUR_COLUMNS, [EUR_DESCRIPTION, EUR_DESCRIPTION2])
    annotate_sas = vcf_annotate(annotate_eur, SAS_AC, SAS_COLUMNS, [SAS_DESCRIPTION, SAS_DESCRIPTION2])
    annotate_afr = vcf_annotate(annotate_sas, AFR_AC, AFR_COLUMNS, [AFR_DESCRIPTION, AFR_DESCRIPTION2])
    annotate_amr = vcf_annotate(annotate_afr, AMR_AC, AMR_COLUMNS, [AMR_DESCRIPTION, AMR_DESCRIPTION2])
    annotate_eas = vcf_annotate(annotate_amr, EAS_AC, EAS_COLUMNS, [EAS_DESCRIPTION, EAS_DESCRIPTION2])

    with open(output_file, 'w') as f:
        f.write(annotate_eas.communicate()[0])


def annotate_hgmd(input_file, output_file):
    """
    Annotate VCF file with known HGMD mutation, gene, and site information.

    Args:
        input_file (str): input VCF file path
        output_file (str): output VCF file path (HGMD information added to file)

    """

    input_text = subprocess.Popen(['cat', input_file], stdout=subprocess.PIPE)

    annotate_hgmd_sites = vcf_annotate(input_text, HGMD, HGMD_SITES_COLUMNS, [HGMD_SITES_DESCRIPTION])
    annotate_hgmd_mut = vcf_annotate(annotate_hgmd_sites, HGMD, HGMD_MUTATIONS_COLUMNS, [HGMD_MUTATIONS_DESCRIPTION])
    annotate_hgmd_genes = vcf_annotate(annotate_hgmd_mut, HGMD_GENES, HGMD_GENES_COLUMNS, [HGMD_GENES_DESCRIPTION])

    with open(output_file, 'w') as f:
        f.write(annotate_hgmd_genes.communicate()[0])


def vcf_annotate(input_text, annotation_source, columns, descriptions):
    """
    Annotate file text using VCF annotate. Input and outputs are subprocess Popen objects for piping.

    Args:
        input_text (Pipe): Popen object to pipe input text
        annotation_source (str): path to VCF annotate formatted txt file with reference variants for annotation
        columns (str): columns contained in annotation_source file
        descriptions (list of str): list of description text inputs for VCF annotate.

    Returns:
        Pipe: object (subprocess) to output pipe of this annotation

    """

    vcf_annotate_path = '/home/unix/mlek/vcftools_0.1.9/bin/vcf-annotate'

    commands = [vcf_annotate_path,
                '-a', annotation_source,
                '-c', columns]

    for description in descriptions:
        commands = commands + ['-d', description]

    pipe = subprocess.Popen(commands, stdin=input_text.stdout, stdout=subprocess.PIPE)
    return pipe


def annotate_dbsnp(input_file, output_file):
    """
    Annotate VCF file with known DBSNP IDs.

    Args:
        input_file (str): input VCF file path
        output_file (str): output VCF file path (DBSNP information added to file)

    """

    input_text = subprocess.Popen(['cat', input_file], stdout=subprocess.PIPE)

    annotate_dbsnp = vcf_annotate(input_text, DBSNP, DBSNP_COLUMNS, [DBSNP_DESCRIPTION])

    with open(output_file, 'w') as f:
        f.write(annotate_dbsnp.communicate()[0])


def annotate_snpeff(input_file, output_file):
    """
    Annotate VCF file using SNPEff.

    Args:
        input_file: input VCF file path
        output_file: output VCF file path (SNPEff annotation added to file)

    """

    snpeff_path = '/seq/software/picard/current/bin/snpEff.jar'
    config_path = '/seq/references/Homo_sapiens_assembly19/v1/snpEff/Homo_sapiens_assembly19.snpEff.config'

    pipe = subprocess.Popen(['java',
                            '-XX:GCTimeLimit=50',
                            '-XX:GCHeapFreeLimit=10',
                            '-Xmx4000m',
                            '-jar', snpeff_path,
                            '-v',
                            '-onlyCoding', 'true',
                            '-c', config_path,
                            '-i', 'vcf',
                            '-o', 'vcf',
                            'GRCh37.64',
                            input_file], stdout=subprocess.PIPE)

    output_text = pipe.communicate()[0]

    with open(output_file, 'w') as f:
        f.write(output_text)