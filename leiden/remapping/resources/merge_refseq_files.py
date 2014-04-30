from leiden_cleanup.input_output import file_io

with open('genes.refGene', 'r') as infile1, open('corrected_genes.refGene', 'r') as infile2:
    file_text1 = infile1.read().splitlines()
    file_text2 = infile2.read().splitlines()

    print len(file_text2)
    print len(file_text1)
    for item in file_text1:
        file_text2.append(item)
    print len(file_text2)
    my_text = set(file_text2)

new_file = []
for lines in my_text:
    new_file.append(lines.split('\t'))

print len(new_file)
file_io.write_table_to_file('combined_genes.refGene', new_file)


