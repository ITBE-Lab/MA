
from setuptools import setup, Extension
 
setup(
    name='Aligner',
    version='0.0.1',
    description='An aligner for Nucleotide Sequences.',
    author='Markus Schmidt',
 
    author_email='markus.rainer.schmidt@gmail.com',
 
    url='http://itbe.hanyang.ac.kr/',
    include_package_data=True,
 
    zip_safe=False,
 
    ext_modules=[
        Extension('Aligner',
                 ['src/yourmodule/foo.cpp',
                   'src/yourmodule/io/someio.cpp',
                   'src/yourmodule/bar.cpp'],
                  include_dirs=['inc'], # assuming your project include files are there
                  #library_dirs=['/somelib/dir'], # optional
                  libraries=['boost_python'], # those are the linked libs
                  extra_compile_args=['-g'] # some other compile args
                  ),
        ]
)

