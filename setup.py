
from setuptools import setup, Extension
 
setup(
    name='LAuS',
    version='0.0.2',
    description='An aligner for Nucleotide Sequences.',
    author='Markus Schmidt',
 
    author_email='markus.rainer.schmidt@gmail.com',
 
    url='http://itbe.hanyang.ac.kr/',
    include_package_data=True,
 
    zip_safe=False,
 
    #license='MIT', #TODO: lookup correct license

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
        #'License :: OSI Approved :: MIT License', #TODO: lookup correct license


        # Specify the Python versions you support here. In particular, ensure 
        # that you indicate whether you support Python 2, Python 3 or both. 
        'Programming Language :: Python :: 2.7', 
    ], 

    keywords='aligner nucleotide sequence',

    ext_modules=[
        Extension('LAuS',
            ['src/aligner.cpp', #TODO: check if this can be automated to pick all .cpp files
            'src/configure.cpp',
            'src/container.cpp',
            'src/FM_index.cpp',
            'src/intervalTree.cpp',
            'src/ksw.c',
            'src/module.cpp',
            'src/nucSeq.cpp',
            'src/pack.cpp',
            'src/segmentation.cpp',
            'src/sequence.cpp',
            'src/support.cpp',
            'src/threadPool.cpp',],
            include_dirs=['inc', 'usr/include', '/opt/dev/boost_1_60_0'], # assuming your project include files are there
            library_dirs=[
                '/usr/lib',
                '/opt/dev/boost_1_60_0/stage/lib', 
                '/opt/dev/lib',
            ], # optional
            libraries=[
                'boost_python', 
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
                'boost_iostreams',
                ], # those are the linked libs
            extra_compile_args=[
                '-std=c++11'
                ] # some other compile args
            ),
        ]
)