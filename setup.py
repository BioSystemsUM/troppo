from setuptools import setup, find_packages

setup(
    name = 'TROPPO',
    version = '0.0.1',
    package_dir = {'':'src'},
    packages = find_packages('src'),
    install_requires = [#"numpy",
                        #scipy",
                        #"pandas",
						"cobamp"],

    author = 'Jorge Ferreira & Vítor Vieira',
    author_email = 'jorge.ferreira@ceb.uminho.pt',
    description = 'TROPPO - Tissue-specific RecOnstruction and Phenotype Prediction using Omics data',
    license = 'GNU General Public License v3.0',
    keywords = 'pathway analysis metabolic model',
    url = 'https://github.com/BioSystemsUM/troppo',
    long_description = open('README.rst').read(),
    classifiers = [
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3.5',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Topic :: Software Development :: Libraries :: Python Modules'
    ],
)