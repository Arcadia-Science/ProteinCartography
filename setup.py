from setuptools import setup

setup(
    name='ProteinCartography',
    url='https://github.com/Arcadia-Science/gene-family-cartography',
    author='Dennis Sun',
    author_email='dennis.sun@arcadiascience.com',
    packages=['ProteinCartography'],
    install_requires=['pandas', 'biopython', 'snakemake', 'matplotlib', 
                      'plotly', 'umap-learn', 'leidenalg', 'scanpy', 'scikit-learn'],
    version='0.0.1',
    license='MIT',
    description='Builds maps of protein space from structures.',
    scripts=['ProteinCartography/aggregate_lists.py',
             'ProteinCartography/esmfold_apiquery.py',
             'ProteinCartography/extract_foldseekhits.py',
             'ProteinCartography/fetch_accession.py',
             'ProteinCartography/foldseek_apiquery.py',
             'ProteinCartography/foldseek_clustering.py',
             'ProteinCartography/make_dummies.py',
             'ProteinCartography/map_refseqids.py',
             'ProteinCartography/run_blast.py']
    # We will also need a readme eventually (there will be a warning)
    # long_description=open('README.txt').read(),
)
