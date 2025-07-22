from setuptools import setup
from Cython.Build import cythonize
from setuptools.extension import Extension

extensions = [
    Extension(
        "gauss_avx2",
        ["gauss_avx2.pyx", "avx2_core.c"],
        extra_compile_args=["-O3", "-mavx2", "-march=native"],
        extra_link_args=["-mavx2"],
    )
]

setup(
    ext_modules=cythonize(extensions)
)