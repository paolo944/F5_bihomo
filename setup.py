from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy

ext = Extension(
    "gauss_avx2",
    sources=["gauss_avx2.pyx", "avx2_core.c"],
    include_dirs=[numpy.get_include()],
    extra_compile_args=["-O3", "-mavx2", "-march=native"],
    extra_link_args=["-mavx2"],
    language="c"
)

setup(
    name="gauss_avx2",
    ext_modules=cythonize([ext]),
)
