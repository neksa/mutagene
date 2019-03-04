import setuptools
from mutagene.version import __version__

setuptools.setup(name='mutagene',
                 version=__version__,
                 description='Mutational analysis with Python',
                 long_description=open('README.md').read().strip(),
                 author='Alexander Goncearenco',
                 author_email='alexandr.goncearenco@nih.gov',
                 url='http://www.ncbi.nlm.nih.gov/projects/mutagene/',
                 py_modules=['mutagene'],
                 install_requires=[],
                 test_suite='nose.collector',
                 tests_require=['nose'],
                 scripts=['bin/mutagene'],
                 include_package_data=True,
                 license='MIT License',
                 zip_safe=False,
                 keywords='bioinformatics cancer mutations',
                 classifiers=[
                        'Development Status :: 4 - Beta',
                        'License :: Public Domain',
                        'Intended Audience :: Science/Research',
                        'License :: OSI Approved :: MIT License',
                        'Programming Language :: Python :: 3 :: Only',
                        'Topic :: Scientific/Engineering :: Bio-Informatics',
                        'Topic :: Scientific/Engineering :: Medical Science Apps.',
                        'Environment :: Console'
                 ])
