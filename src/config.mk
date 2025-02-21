# Choose your compiler by setting FC to either 'ifort' or 'gfortran'
FC = gfortran

# Compiler-specific flags
ifeq ($(FC),ifort)
    FFLAGS_DEBUG = -g
    FFLAGS_OPT = -O2
    LDFLAGS = -g
else ifeq ($(FC),gfortran)
    FFLAGS_DEBUG = -g -fdefault-real-8 -fdefault-double-8
    FFLAGS_OPT = -O2 -fdefault-real-8 -fdefault-double-8
    LDFLAGS = -g
else
    $(error Unsupported compiler: $(FC). Use either 'ifort' or 'gfortran')
endif