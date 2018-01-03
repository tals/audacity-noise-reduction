from distutils.core import setup, Extension
from Cython.Build import cythonize
from glob import glob

sources = [x for x in glob('noisereduction/*.cpp') if not x.endswith('main.cpp')]
setup(ext_modules = cythonize(Extension(
           "noisereduction.pyx",
           sources=sources,
           language="c++",
           extra_compile_args=["-std=c++11"],
           extra_link_args=["-std=c++11"]

      )))