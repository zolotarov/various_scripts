import csv

five_LF = set()
clark = set()
glab = set()

with open("epidermis_up.csv", 'r') as gene_csv:
    gene_reader = csv.reader(gene_csv, delimiter=',', quotechar='|')
    gene_reader.next()  # skip the header row
    for row in gene_reader:
        five_LF.add(row[0])
        clark.add(row[1])
        glab.add(row[2])


common = glab & five_LF

difference = glab - (clark | five_LF)  # genes unique to the first set, not found in sets 2 and 3

for item in difference:
    print item

