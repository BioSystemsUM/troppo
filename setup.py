from setuptools import setup, find_packages

setup(
    name='troppo',
    version='0.0.6',
    package_dir={'': 'src'},
    packages=find_packages('src'),
    install_requires=["cobamp==0.2.1"
                      "numpy==1.21.6"
                      "cobra==0.26.2"
                      "pandas==1.3.5"
                      "troppo==0.0.6"
                      "pathos==0.2.9"
                      "scipy==1.7.3"
                      "requests==2.28.2"
                      "matplotlib==3.5.2"
                      "setuptools==65.6.3"],

    author='Jorge Ferreira & Vítor Vieira',
    author_email='jorge.ferreira@ceb.uminho.pt',
    description='TROPPO - Tissue-specific RecOnstruction and Phenotype Prediction using Omics data',
    license='GNU General Public License v3.0',
    keywords='pathway analysis metabolic model',
    url='https://github.com/BioSystemsUM/troppo',
    long_description=open('README.rst').read(),
    classifiers=[
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3.5',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Topic :: Software Development :: Libraries :: Python Modules'
    ],
)
