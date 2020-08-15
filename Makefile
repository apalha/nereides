# make file that compiles all fortran modules for nereides
# using f2py
all: _mesher _math_tools

_mesher:
	f2py -c --fcompiler=gnu95 -m _mesher _mesher.f90 -DF2PY_REPORT_ON_ARRAY_COPY=1

_math_tools:
	f2py -c --fcompiler=gnu95 -m _math_tools _math_tools.f90 -DF2PY_REPORT_ON_ARRAY_COPY=1
