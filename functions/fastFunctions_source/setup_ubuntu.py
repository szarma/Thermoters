from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np

ext_modules=[
    Extension("fastFunctions",
              ["fastFunctions.pyx"],
              libraries   = ["m"],
              include_dirs = [np.get_include()],
#               extra_compile_args = [
#               "-O3",
#               "-ffast-math",
#               "-march=native"
#               ],
#               extra_link_args=['-fopenmp']
              ) 
]
# 
setup( 
  name = "cythonized functions for promoter prediction",
  cmdclass = {"build_ext": build_ext},
  include_dirs = [np.get_include()],
  ext_modules = ext_modules
)

# include_dir=[np.get_include()] to Extension(..., include_dir=[np.get_include()])
# in addition to setup(...,include_dir=[np.get_include()])

