from setuptools import find_packages
from numpy.distutils.core import setup#, Extension
#from pathlib import Path

#BASE_DIR = Path(__file__).resolve().parent

#with open('requirements.txt') as f:
#    install_requires = f.read().splitlines()

#ext1 = Extension(name='m3l.dr_module',
#                 sources=['src/m3l/dr_module.f90'],
#                 f2py_options=['--quiet'],
#                 )

if __name__ == "__main__":

    setup(
    #    install_requires=install_requires,
        name = "m3l",
        version = "2023.05.dev0",
        url = "https://github.com/flavianowilliams/M3L",
        package_dir={"":"src"},
        packages=find_packages(where="src"),
        description = "Machine Learning and Molecular Modelling",
        author = "Flaviano Williams Fernandes",
        author_email = "flaviano.fernandes@ifpr.edu.br",
        ext_modules=[]
    )
