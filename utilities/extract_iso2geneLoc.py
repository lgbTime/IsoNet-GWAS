import re
import sys
combined_gtf = sys.argv[1]
with open(combined_gtf, 'r') as gff:
    lines = gff.readlines()
for line in lines:
    iso = line.strip().split("transcript_id ")[1].split(";")[0].replace('"', '')
    geneLOC = re.search(r"LOC\d+", line.strip())
    geneOther = re.search(r"gene_name \w+;",line.strip())
    chrome = line.strip().split("\t")[0]
    starts = line.strip().split("\t")[3]
    ends =   line.strip().split("\t")[4]
    if geneLOC:
        print(chrome,starts,ends,iso,geneLOC[0],sep="\t")
    elif geneOther:
        print(chrome,starts,ends,iso,geneOther.split(";")[0],sep="\t")
    else:
        print(chrome,starts,ends,iso,iso,sep="\t")


