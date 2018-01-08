
from setuptools import setup, Extension
from os import listdir

src_folders = list(map(lambda x: "src/"+x, listdir("src")))

src_files = []
for folder in src_folders:
    src_files.extend(list(map(lambda x: folder+"/"+x, listdir(folder))))

setup(
    name='MABS',
    version='0.0.4',
    description='An aligner for Nucleotide Sequences.',
    author='Markus Schmidt',
 
    author_email='markus.rainer.schmidt@gmail.com',
 
    url='http://itbe.hanyang.ac.kr/',
    include_package_data=True,
 
    zip_safe=False,
    python_requires='>=3',
 
    license='MIT', #@TODO: lookup correct license

    classifiers=[ 
        # How mature is this project? Common values are 
        #   3 - Alpha 
        #   4 - Beta 
        #   5 - Production/Stable 
        'Development Status :: 3 - Alpha', 


        # Indicate who your project is intended for 
        'Intended Audience :: Developers', 
        'Topic :: Scientific/Engineering :: Bio-Informatics', 

        # Pick your license as you wish (should match "license" above) 
        'License :: OSI Approved :: MIT License', #@TODO: lookup correct license


        # Specify the Python versions you support here. In particular, ensure 
        # that you indicate whether you support Python 2, Python 3 or both. 
        'Programming Language :: Python :: 3.5', 
    ], 

    keywords='aligner nucleotide sequence',

    ext_modules=[
        Extension(
            'libMABS',
            src_files,
            include_dirs=['inc', '/opt/dev/boost_1_65_1'], # project include files
            library_dirs=[
                '/usr/lib',
                '/opt/dev/boost_1_65_1/stage/lib', 
                '/opt/dev/lib',
            ], # optional
            libraries=[
                'boost_python3-mt', 
                'dl',
                'rt',
                'z',
                'boost_system-mt',
                'boost_thread-mt',
                'boost_log-mt',
                'boost_log_setup-mt' ,
                'boost_filesystem-mt' ,
                'boost_program_options-mt' ,
                'boost_regex-mt',
                'boost_iostreams-mt',
                ], # those are the linked libs
            extra_compile_args=[
                '-std=c++11',
                '-DBOOST_ALL_DYN_LINK'
                ] # some other compile args
            ),
        ],
    packages=["MABS"]
)