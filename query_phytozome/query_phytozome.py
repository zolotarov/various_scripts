# -*- coding: UTF-8 -*-

"""
Original script: https://gist.github.com/384470#file-gistfile1-py

This script queries the phytozome database
In the current form, you can either obtain 1000 bp of the promoter sequences
or peptide sequences

To modify the script with other queries, go to http://www.phytozome.net/biomart/martview
create a custom query and then click on XML, to get an example code
"""


import urllib as U
from collections import namedtuple
Species_info = namedtuple('Species_info', 'species_name shorthand value')

list_of_species = [Species_info('Aquilegia coerulea', 'aquco', '195'), Species_info('Arabidopsis lyrata', 'araly', '107'), Species_info('Arabidopsis thaliana', 'arath', '167'), Species_info('Brachypodium distachyon', 'bradi', '192'), Species_info('Brassica rapa', 'brara', '197'), Species_info('Capsella rubella ', 'capru', '183'), Species_info('Carica papaya', 'carpa', '113'), Species_info('Chlamydomonas reinhardtii', 'chlre', '236'), Species_info('Citrus clementina', 'citcl', '182'), Species_info('Citrus sinensis', 'citsi', '154'), Species_info('Coccomyxa subellipsoidea C-169', 'cocsu', '227'), Species_info('Cucumis sativus', 'cucsa', '122'), Species_info('Eucalyptus grandis', 'eucgr', '201'), Species_info('Fragaria vesca', 'frave', '226'), Species_info('Glycine max', 'glyma', '189'), Species_info('Gossypium raimondii', 'gosra', '221'), Species_info('Linum usitatissimum', 'linus', '200'), Species_info('Malus domestica', 'maldo', '196'), Species_info('Manihot esculenta', 'manes', '147'), Species_info('Medicago truncatula', 'medtr', '198'), Species_info('Micromonas pusilla CCMP1545', 'micpu', '228'), Species_info('Micromonas pusilla RCC299', 'micpv', '229'), Species_info('Mimulus guttatus', 'mimgu', '140'), Species_info('Oryza sativa', 'orysa', '204'), Species_info('Ostreococcus lucimarinus', 'ostlu', '231'), Species_info('Panicum virgatum', 'panvi', '202'), Species_info('Phaseolus vulgaris', 'phavu', '218'), Species_info('Physcomitrella patens', 'phypa', '152'), Species_info('Populus trichocarpa', 'poptr', '210'), Species_info('Prunus persica', 'prupe', '139'), Species_info('Ricinius communis', 'ricco', '119'), Species_info('Selaginella moellendorffii', 'selmo', '91'), Species_info('Setaria italica', 'setit', '164'), Species_info('Solanum lycopersicum', 'solly', '225'), Species_info('Solanum tuberosum', 'soltu', '206'), Species_info('Sorghum bicolor', 'sorbi', '79'), Species_info('Thellungiella halophila', 'theha', '173'), Species_info('Theobroma cacao', 'theca', '233'), Species_info('Vitis vinifera', 'vitvi', '145'), Species_info('Volvox carteri', 'volca', '199'), Species_info('Zea mays', 'zeama', '181')]


xml_promoter = """
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "zome_mart" formatter = "FASTA" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >

	<Dataset name = "phytozome" interface = "default" >
		<Filter name = "upstream_flank" value = "1000"/>
		<Filter name = "organism_id" value = "%s"/>
		<Attribute name = "gene_name1" />
		<Attribute name = "transcript_name1" />
		<Attribute name = "coding_gene_flank" />
		<Attribute name = "organism_name" />
	</Dataset>
</Query>
""".replace("\n", "")

xml_peptide = """
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "zome_mart" formatter = "FASTA" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
			
	<Dataset name = "phytozome" interface = "default" >
		<Filter name = "organism_id" value = "%s"/>
		<Attribute name = "peptide_sequence" />
		<Attribute name = "gene_name1" />
		<Attribute name = "transcript_name1" />
		<Attribute name = "organism_name" />
	</Dataset>
</Query>
""".replace("\n", "")

url = "http://www.phytozome.net/biomart/martservice"
# check available filters and attributes at following address
filters = "%s?type=filters;dataset=phytozome" % url
attributes = "%s?type=attributes;dataset=phytozome" % url


class PhytozomeMart(object):

    def __init__(self, xml, url, filters):
        # xml = xml % ",".join(str(x) for x in filter)
        xml = xml % filters
        self.xml = xml
        self.url = url

    def send_query(self):
        data = U.urlencode({"query": self.xml})
        return U.urlopen(self.url + "?" + data).read()


# if __name__ == '__main__':
# 	mart_output = open('test.fas', 'w')
# 	mart = PhytozomeMart(xml, url, filters = 195)
# 	mart_output.write(mart.send_query())

# for item in list_of_species:
# 	mart_output = open('%s_prom_1K.fas' % item.shorthand, 'w')
# 	print "Downloading promoters for %s" % item.species_name
# 	mart = PhytozomeMart(xml, url, filters = item.value)
# 	mart_output.write(mart.send_query())

# for item in list_of_species:
# 	mart_output = open('%s_pep.fas' % item.shorthand, 'w')
# 	mart = PhytozomeMart(xml_peptide, url, item.value)
# 	mart_output.write(mart.send_query())

def download_sequences(kind):
	if kind == "promoters":
		for item in list_of_species:
			mart_output = open('%s_prom_1K.fas' % item.shorthand, 'w')
			mart = PhytozomeMart(xml_promoter, url, item.value)
			mart_output.write(mart.send_query())
	elif kind == "peptides":
		for item in list_of_species:
			mart_output = open('%s_pep.fas' % item.shorthand, 'w')
			mart = PhytozomeMart(xml_peptide, url, item.value)
			mart_output.write(mart.send_query())

download_sequences("promoters")

