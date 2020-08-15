# import platform to know the type of processor the Nereides in running on
# 32bit or 64bit
import platform

if platform.architecture()[0] == "32bit":
    ## Type of \c int used for indexing sparse matrices
    # dependent on the processor architecture. int32 is used for 32bit machines.
    index_INT_TYPE = "int32"
else:
    ## Type of \c int used for indexing sparse matrices
    # dependent on the processor architecture. int64 is used for 64bit machines.
    index_INT_TYPE = "int64"