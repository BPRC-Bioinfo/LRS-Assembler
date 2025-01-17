import json
import pandas as pd
import sys

inFile = sys.argv[1]

# Open and read the JSON file
accession = []
chrName = []

with open(inFile) as f:
    line_content = [json.loads(line) for line in f.readlines()]
    for line in line_content:
        if line['role'] == 'assembled-molecule':
            accession.append(line['refseqAccession'])
            # Check if the chromosome name is a digit and less than 10
            if line['chrName'].isdigit() and int(line['chrName']) < 10:
                chrName.append(f'{int(line["chrName"]):02d}')  # Format with leading zero
            else:
                chrName.append(line['chrName'])

df = pd.DataFrame({'accession': accession, 'chromosome': chrName})
print(df.to_string(index=False))

