from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np

compile_args = ["-O3", "-mavx2", "-march=native", "-mtune=native", "-funroll-loops", "-flto", "-fopenmp"]
link_args = ["-fopenmp"]

extensions = [
    Extension(
        "gauss_gf2",
        sources=["gauss_gf2.pyx", "gauss_core.c"],
        include_dirs=[np.get_include(), "."],
        extra_compile_args=compile_args,
        extra_link_args=link_args,
        language="c",
    ),
    Extension(
        "pack_utils",
        sources=["pack_utils.pyx"],  # ✅ correction ici
        include_dirs=[np.get_include(), "."],
        extra_compile_args=compile_args,
        extra_link_args=link_args,
        language="c",
    )
]

setup(
    ext_modules=cythonize(extensions, language_level=3),
)
