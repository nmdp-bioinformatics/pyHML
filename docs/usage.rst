=====
Usage
=====

To use pyHML in a project::

    from pyhml import pyhml
    hml_file = "hml_example.xml"
    hmlparser = pyhml.HmlParser()
    hml = hmlparser.parse(hml_file)
    pandasdf = pyhml.toDF(hml)


