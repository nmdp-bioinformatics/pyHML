=====
Usage
=====

To use pyHML in a project::

    import pyhml 
	hmlparser = pyhml.HmlParser()
    hml = hmlparser.parse("hml_example.xml")
    pandasdf = hml.toPandas()

    # Ouput the HML data as a IPD-IMGT/HLA .dat file for each subject
    hml.tobiotype("output/directory", dtype='imgt', by='subject')

    # Output the whole HML file as one fasta file
    hml.tobiotype("output/directory", dtype='fasta', by='file')

    # Defaults to dtype='fasta' and by='subject'
    hml.tobiotype("output/directory")

