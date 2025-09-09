import json
import pandas as pd
import sys

inFile = sys.argv[1]

accessions = []
chromosomes = []

with open(inFile, 'r') as f:
    for raw in f:
        line = json.loads(raw)
        
        # Prefer refseq, fallback to genbank
        acc = line.get('refseqAccession') or line.get('genbankAccession')
        if not acc:
            continue  # skip entries with neither accession
        
        accessions.append(acc)
        
        # Chromosome formatting
        chr_name = line.get('chrName', '')
        if chr_name.isdigit() and int(chr_name) < 10:
            chromosomes.append(f'{int(chr_name):02d}')
        else:
            chromosomes.append(chr_name)

# Create and print DataFrame
df = pd.DataFrame({
    'accession': accessions,
    'chromosome': chromosomes
})
print(df.to_string(index=False))