#!/usr/bin/python

"""
Setup script for bioitools
"""

try:
	from setuptools import setup, Extension
except:
	from distutils.core import setup, Extension

#mod1 = Extension('cIO', sources=["fastqc_py/cIO.c"])

setup(name = "bioitools",
	version = "0.1",
	author = "Greg Zynda",
	author_email="gzynda@tacc.utexas.edu",
	license="GNU",
	description="Bioinformatics tools",
	install_requires=['numpy','PyWavelets'],
	packages = ["bioitools"],
	scripts = ["bin/bioitools"],
	test_suite="tests")
	#ext_modules=[mod1],
	#package_data={'bioitools/test_data':['test.fa.fai']})
