from leiden_cleanup.input_output import file_io

with open('genes.refGene', 'r') as infile:
    file_text = infile.read().splitlines()

    new_lines = []
    for line in file_text:
        entries = line.split('\t')

        version = entries[-1]
        del entries[-1]

        entries[1] = entries[1] + '.' + version
        new_lines.append(entries)

    file_io.write_table_to_file('corrected_genes.refGene', new_lines)
