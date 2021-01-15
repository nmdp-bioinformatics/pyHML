#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    pyhml pyHML.
#    Copyright (c) 2017 Be The Match operated by National Marrow Donor Program. All Rights Reserved.
#
#    This library is free software; you can redistribute it and/or modify it
#    under the terms of the GNU Lesser General Public License as published
#    by the Free Software Foundation; either version 3 of the License, or (at
#    your option) any later version.
#
#    This library is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; with out even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
#    License for more details.
#
#    You should have received a copy of the GNU Lesser General Public License
#    along with this library;  if not, write to the Free Software Foundation,
#    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA.
#
#    > http://www.fsf.org/licensing/licenses/lgpl.html
#    > http://www.opensource.org/licenses/lgpl-license.php
#


from setuptools import setup

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'biopython==1.71',
    'pandas==0.20.3',
    'numpy==1.14.2',
    'six==1.11.0',
    'xmlschema==0.9.13',
    'xmltodict==0.11.0',
    'pytz==2018.3',
    'sh==1.12.14',
    'python-dateutil==2.7.2'
]

test_requirements = [
    # TODO: put package test requirements here
]

setup(
    name='pyhml',
    version='0.0.6',
    description="Python HML parser",
    long_description=readme + '\n\n' + history,
    author="Mike Halagan",
    author_email='mhalagan@nmdp.org',
    url='https://github.com/nmdp-bioinformatics/pyHML',
    packages=[
        'pyhml',
        'pyhml.models'
    ],
    package_dir={'pyhml':
                 'pyhml'},
    package_data={'pyhml': ['data/*']},
    install_requires=requirements,
    license="LGPL 3.0",
    zip_safe=False,
    keywords='pyhml',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
    ],
    test_suite='tests',
    tests_require=test_requirements
)
