from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np
import os

ext = Extension(
    "gauss_gf2",
    sources=["gauss_gf2.pyx", "gauss_core.c"],
    include_dirs=[np.get_include()],
    extra_compile_args=["-O3", "-mavx2", "-march=native", "-mtune=native", "-funroll-loops", "-flto", "-fopenmp"],
    extra_link_args=['-fopenmp'],
    language="c",
)

setup(
    name="gauss_gf2",
    ext_modules=cythonize([ext], language_level=3),
)