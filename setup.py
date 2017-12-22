import sys, os
from os.path import relpath, join
from setuptools import setup

if sys.argv[-1] == 'setup.py':
    print("To install, run 'python setup.py install'")
    print()

def find_package_data(data_root, package_root):
    files = []
    for root, dirnames, filenames in os.walk(data_root):
        for fn in filenames:
            files.append(relpath(join(root, fn), package_root))
    return files


descr = """

""" #TODO

setup(
    name = 'rdkommtools',
    version = '0.0.1',
    description = 'RDKit OpenMM Tools',
    long_description = descr,
    url = 'https://github.com/hjuinj/rdkommtools',
    license = 'MIT',
    author = 'Shuzhe Wang',
    author_email = 'shuwang.wang@phys.chem.ethz.ch',
    platforms = ['Linux-64', 'Mac OSX-64', 'Unix-64'],
    packages = find_packages()+['tests'],
    include_package_data = True,
    zip_safe = False
)
