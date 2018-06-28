import setuptools
from pymutagene.version import Version


setuptools.setup(name='pymutagene',
                 version=Version('0.0.1').number,
                 description='Mutational analysis with python',
                 long_description=open('README.md').read().strip(),
                 author='Alexander Goncearenco',
                 author_email='alexandr.goncearenco@nih.gov',
                 url='http://www.ncbi.nlm.nih.gov/projects/mutagene/',
                 py_modules=['pymutagene'],
                 install_requires=[],
                 license='MIT License',
                 zip_safe=False,
                 keywords='bioinformatics cancer mutations',
                 classifiers=['Mutations', ])
