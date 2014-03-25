# -*- coding: utf-8 -*- 

import requests

url = "http://www.phytozome.net/biomart/martservice?query=%s"
xml = """
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "zome_mart" formatter = "FASTA" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
			
	<Dataset name = "phytozome" interface = "default" >
		<Filter name = "upstream_flank" value = "1000"/>
		<Filter name = "organism_id" value = "189"/>
		<Attribute name = "gene_name1" />
		<Attribute name = "transcript_name1" />
		<Attribute name = "gene_flank" />
	</Dataset>
</Query>""".replace("\n", "")

# print xml

# query = {'query': xml}
headers = {'Content-Type': 'application/xml'}
r = requests.get(url % xml)
print r.url
# print r.status_code
# print r.headers
# print r.content
# print r.text
with open('temp.fas', 'w') as file:
	file.write(r.text)

