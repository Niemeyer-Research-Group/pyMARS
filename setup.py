"""Setup script to build module.
"""
#!/usr/bin/env python

from setuptools import setup
from codecs import open
from os import path
import sys

here = path.abspath(path.dirname(__file__))

version = {}
with open(path.join(here, 'pymars', '_version.py')) as version_file:
    exec(version_file.read(), version)
__version__ = version['__version__']

with open(path.join(here, 'README.md')) as readme_file:
    readme = readme_file.read()

with open(path.join(here, 'CHANGELOG.md')) as changelog_file:
    changelog = changelog_file.read()

with open(path.join(here, 'CITATION.md')) as citation_file:
    citation = citation_file.read()

long_description = readme + '\n\n' + changelog + '\n\n' + citation

install_requires = [
    'argparse',
    'h5py',
    'numpy',
    'cantera>=2.3.0',
    'networkx',
    'pytables'
]

tests_require = [
    'pytest',
    'pytest-cov',
]

needs_pytest = {'pytest', 'test', 'ptr'}.intersection(sys.argv)
setup_requires = ['pytest-runner'] if needs_pytest else []

setup(
    name='pymars',
    version=__version__,
    description=(
        'Python-based chemical kinetic Model Automatic '
        'Reduction Software (pyMARS)'
        ),
    long_description=long_description,
    url='https://github.com/Niemeyer-Research-Group/pyMARS',
    author='Phillip Mestas, Parker Clayton, Kyle Niemeyer',
    author_email='kyle.niemeyer@gmail.com',
    license='MIT',
    keywords=['chemical kinetics'],
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
    ],
    packages=['pymars', 'pymars.tests'],
    package_dir={'pymars': 'pymars'},
    include_package_data=True,
    entry_points={
        'console_scripts': ['pymars=pymars.__main__:main']
        },
    install_requires=install_requires,
    tests_require=tests_require,
    setup_requires=setup_requires,
    python_requires='~=3.5',
    zip_safe=False,
    )
