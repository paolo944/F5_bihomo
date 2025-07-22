from setuptools import setup, Extension
from Cython.Build import cythonize

ext = Extension(
    name="gauss_avx2",
    sources=["gauss_avx2.pyx", "avx2_core.c"],
    extra_compile_args=["-mavx2", "-O3"],
)

setup(
    ext_modules=cythonize([ext])
)