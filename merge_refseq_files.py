import lovd.database.utilities

with open('lovd/remapping/resources/genes.refGene', 'r') as infile1, open('lovd/remapping/resources/original.refGene', 'r') as infile2:
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
lovd.database.utilities.write_table_to_file('combined_genes.refGene', new_file)


