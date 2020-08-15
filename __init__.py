# expose definitions from modules to this package

# import platform to know the type of processor the Nereides in running on
# 32bit or 64bit
import platform

##if platform.architecture()[0] == "32bit":
##    ## Type of \c int used for indexing sparse matrices
##    # dependent on the processor architecture. int32 is used for 32bit machines.
##    index_INT_TYPE = "int32"
##else:
##    ## Type of \c int used for indexing sparse matrices
##    # dependent on the processor architecture. int64 is used for 64bit machines.
##    index_INT_TYPE = "int64"

# import mesh modules
import mesher
import gauss_quadrature
import math_tools
import fictitious_domain
import settings
import physics

# import exceptions
import NereidesExceptions