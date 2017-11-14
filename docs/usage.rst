=====
Usage
=====

To use pyHML in a project::

    import pyhml 
	hmlparser = pyhml.HmlParser()
    hml = hmlparser.parse("hml_example.xml")
    pandasdf = pyhml.toDF(hml)


