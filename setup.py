import setuptools

setuptools.setup(
    name='wgbs_tools',
    version='0.1',
    description='Toolkit to manipulate and analyze Whole Genome Bisulfite Sequencing data',
    packages=setuptools.find_packages(),
    include_package_data=True,
    py_modules=['wgbs_tools'],
    install_requires=[
        'Click',
        'pysam',
        'pyyaml',
    ],
    entry_points='''
        [console_scripts]
        wgbs_tools=wgbs_tools_cli:cli
    '''
)