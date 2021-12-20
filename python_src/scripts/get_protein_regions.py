import pandas as pd
import sys

def get_gene_product(x):
    """
    parse x.info and retrieve the gene (if present)
    and product (ie protein)
    """
    gene = None
    product = None
    for item in x.info.split(';'):
        key, value = item.split('=')
        if key == 'gene':
            gene = value
        elif key == 'product':
            product = value
    if not gene:
        gene = 'ORF1ab'
    x['gene'], x['product'] = gene, product
    return x


annotations = sys.argv[1]
output = sys.argv[2]

# 0 - contig
# 2 - type of region
# 3 - start
# 4 - end
# 8 - info string (; delimited)
data = pd.read_csv(annotations, sep='\t', compression='gzip',
                   comment='#', usecols=[0,2,3,4,8],
                   names=['contig', 'region_type', 'start', 'end', 'info'])
data = data[data.region_type.str.contains('CDS') &
            ~data["info"].str.contains('gene=ORF1ab')]

data = data.apply(
    get_gene_product,
    axis=1,
    result_type='expand')

data[['contig', 'start', 'end', 'gene', 'product']] \
    .to_csv(output, sep='\t', index=False, header=False)
