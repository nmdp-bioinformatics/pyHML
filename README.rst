===============================
pyHML
===============================


.. image:: https://img.shields.io/pypi/v/pyhml.svg
        :target: https://pypi.python.org/pypi/pyhml

.. image:: https://img.shields.io/travis/mhalagan-nmdp/pyhml.svg
        :target: https://travis-ci.org/mhalagan-nmdp/pyhml

.. image:: https://readthedocs.org/projects/pyhml/badge/?version=latest
        :target: https://pyhml.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. image:: https://pyup.io/repos/github/mhalagan-nmdp/pyhml/shield.svg
     :target: https://pyup.io/repos/github/mhalagan-nmdp/pyhml/
     :alt: Updates


Python HML parser


* Free software: LGPL 3.0
* Documentation: https://pyhml.readthedocs.io.


Features
--------

.. code-block:: python3

    from pyhml.pyhml import HmlParser, toDF, tobiotype
    hml_file = "hml_example.xml"
    hmlparser = HmlParser()
    hml = hmlparser.parse(hml_file)
    outdir = 'output/directory'

    # Print out each subject in fasta format
    tobiotype(hml, outdir, dtype='fasta', by='subject')

    # Print out the full HML file in IMGT dat file format
    tobiotype(hml, outdir, dtype='imgt', by='hml')

    # Get pandas DF from HML object
    pandasdf = toDF(hml)
    print(pandasdf)

             ID     Locus                             glstring dbversion  \
	0   1367-7150-8     HLA-A        HLA-A*01:01:01+HLA-A*24:02:01    3.14.0   
	1   1367-7150-8     HLA-A        HLA-A*01:01:01+HLA-A*24:02:01    3.14.0   
	2   1367-7150-8     HLA-A        HLA-A*01:01:01+HLA-A*24:02:01    3.14.0   
	3   1367-7150-8     HLA-A        HLA-A*01:01:01+HLA-A*24:02:01    3.14.0   
	4   1367-7150-8     HLA-B        HLA-B*08:01:01+HLA-B*57:01:01    3.14.0   
	5   1367-7150-8     HLA-B        HLA-B*08:01:01+HLA-B*57:01:01    3.14.0   
	6   1367-7150-8     HLA-B        HLA-B*08:01:01+HLA-B*57:01:01    3.14.0   
	7   1367-7150-8     HLA-B        HLA-B*08:01:01+HLA-B*57:01:01    3.14.0   
	8   1367-7150-8     HLA-C        HLA-C*06:02:01+HLA-C*07:01:01    3.14.0   
	9   1367-7150-8     HLA-C        HLA-C*06:02:01+HLA-C*07:01:01    3.14.0   
	10  1367-7150-8     HLA-C        HLA-C*06:02:01+HLA-C*07:01:01    3.14.0   
	11  1367-7150-8     HLA-C        HLA-C*06:02:01+HLA-C*07:01:01    3.14.0   
	12  1367-7150-8  HLA-DPB1  HLA-DPB1*02:01:02+HLA-DPB1*04:01:01    3.14.0   
	13  1367-7150-8  HLA-DPB1  HLA-DPB1*02:01:02+HLA-DPB1*04:01:01    3.14.0   
	14  1367-7150-8  HLA-DRB1  HLA-DRB1*03:01:01+HLA-DRB1*07:01:01    3.15.0   
	15  1367-7150-8  HLA-DRB1  HLA-DRB1*03:01:01+HLA-DRB1*07:01:01    3.15.0   

	                                             sequence  
	0   TTCCTGGATACTCACGACGCGGACCCAGTTCTCACTCCCATTGGGT...  
	1   TTCCCGTCAGACCCCCCCAAGACACATATGACCCACCACCCCATCT...  
	2   TTCCTGGATACTCACGACGCGGACCCAGTTCTCACTCCCATTGGGT...  
	3   GTGCCTGTGTCCAGGCTGGTGTCTGGGTTCTGTGCTCTCTTCCCCA...  
	4   CCATGGTGAGTTTCCCTGTACAAGAGTCCAAGGGGAGAGGTAAGTG...  
	5   GGCCTCTGCGGAGAGGAGCGAGGGGCCCGCCCGGCGAGGGCGCAGG...  
	6   CCATGGTGAGTTTCCCTGTACAAGAGTCCAAGGGGAGAGGTAAGTG...  
	7   GGCCTCTGCGGAGAGGAGCGAGGGGCCCGCCCGGCGAGGGCGCAGG...  
	8   AGGGATCAGGACGAAGTCCCAGGTCCCGGACGGGGCTCTCAGGGTC...  
	9   CGCATCCCCACTTCCCACTCCCATTGGGTGTCGGATATCTAGAGAA...  
	10  AGGGATCAGGACGAAGTCCCAGGTCCCGGACGGGGCTCTCAGGGTC...  
	11  CGCATCCCCACTTCCCACTCCCATTGGGTGTCGGATATCTAGAGAA...  
	12  CCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATT...  
	13  CCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATT...  
	14  CATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA...  
	15  CATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA... 


Install
--------

.. code-block:: bash

    pip install -e 'git+https://github.com/mhalagan-nmdp/pyHML.git#egg=0.0.1'


Credits
---------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage

