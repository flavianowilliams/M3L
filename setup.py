from setuptools import find_packages
from numpy.distutils.core import setup, Extension

#with open('requirements.txt') as f:
#    install_requires = f.read().splitlines()

ext1 = Extension(name='m3l.dr_module',
                 sources=['src/m3l/dr_module.f95'],
                 f2py_options=['--quiet'],
                 )

setup(
#    install_requires=install_requires,
    package_dir={"":"src"},
    packages=find_packages(where="src"),
    ext_modules=[ext1]
)
