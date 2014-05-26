from distutils.core import setup

setup(
    name='BamM',
    version='0.0.1',
    author='Michael Imelfort',
    author_email='mike@mikeimelfort.com',
    packages=['bamm'],
    scripts=['bin/BamM'],
    url='http://pypi.python.org/pypi/BamM/',
    license='GPLv3',
    description='BamM',
    long_description=open('README.md').read(),
    install_requires=[],
)

