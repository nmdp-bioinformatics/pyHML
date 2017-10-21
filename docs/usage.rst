=====
Usage
=====

To use pyHML in a project::

    from pyhml.pyhml import HmlParser, toDF, tobiotype
    hml_file = "hml_example.xml"
    hmlparser = HmlParser()
    hml = hmlparser.parse(hml_file)
   	pandasdf = toDF(hml)


