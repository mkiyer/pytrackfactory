'''
Created on Mar 3, 2011

@author: mkiyer
'''
from distutils.core import setup
from distutils.sysconfig import get_python_inc, get_python_lib
from distutils.extension import Extension
from Cython.Distutils import build_ext
import os

numpy_inc = os.path.join(get_python_lib(), 'numpy/core/include')

#ext_modules = [Extension("pytrackfactory.cintervaltree", ["pytrackfactory/cintervaltree.pyx"], include_dirs=[numpy_inc])]
#               Extension("trackfactory.io.cintervals", ["trackfactory/io/cinterval.pyx"], include_dirs=[numpy_inc]),

ext_modules = [Extension("trackfactory.io.cwiggle", ["trackfactory/io/cwiggle.pyx"]),
               Extension("trackfactory.io.cbedgraph", ["trackfactory/io/cbedgraph.pyx"], include_dirs=[numpy_inc])]

setup(name='pytrackfactory',
      cmdclass={'build_ext': build_ext},
      ext_modules=ext_modules)
