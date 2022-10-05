from setuptools import find_packages
from numpy.distutils.core import setup, Extension

ext1 = Extension(name='m3l.dr_module',
                 sources=['src/m3l/dr_module.f95'],
                 f2py_options=['--quiet'],
                 )

setup(
    package_dir={"":"src"},
    packages=find_packages(where="src"),
    ext_modules=[ext1]
)
