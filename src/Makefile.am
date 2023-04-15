
if BUNDLED_ISL
    ISL_LA = $(top_builddir)/isl/libisl.la
endif
if BUILD_ISL
    ISL_LA = $(ISL_BUILDDIR)/libisl.la
endif

if BUNDLED_BARVINOK
    BARVINOK_LA = $(top_builddir)/barvinok/libbarvinok.la
endif
if BUILD_BARVINOK
    BARVINOK_LA = $(BARVINOK_BUILDDIR)/libbarvinok.la
endif

if BUNDLED_PET
    PET_LA = $(top_builddir)/pet/libpet.la
endif
if BUILD_PET
    PET_LA = $(PET_BUILDDIR)/libpet.la
endif

if DEBUG
    AM_CFLAGS = -g3 -O0
    AM_CXXFLAGS = -g3 -O0
else
    AM_CFLAGS = -O3 -DNDEBUG
    AM_CXXFLAGS = -O3 -DNDEBUG
endif

FORCE:
../isl/libisl.la: FORCE
	cd ../isl; $(MAKE) $(AM_MAKEFLAGS) libisl.la
../barvinok/libbarvinok.la: FORCE
	cd ../barvinok; $(MAKE) $(AM_MAKEFLAGS) libbarvinok.la
../pet/libpet.la: FORCE
	cd ../pet; $(MAKE) $(AM_MAKEFLAGS) libpet.la

AM_CPPFLAGS = @ISL_CFLAGS@ @BARVINOK_CFLAGS@ @PET_CFLAGS@ 
AM_LDFLAGS =  -lcplex -lconcert -lilocplex -lm -lpthread -ldl -lstdc++ 
LIB_ISL = $(ISL_LA) @ISL_LIBS@
LIB_BARVINOK = $(BARVINOK_LA) @BARVINOK_LIBS@
LIB_PET = $(PET_LA) @PET_LIBS@
LIB_CPLEX = $(CPLEX_LA) @CPLEX_LIBS@
LINK = $(CXXLINK)

# build the bullseye library
noinst_LTLIBRARIES = libbullseye.la
libbullseye_la_SOURCES = \
    Program.cpp \
    BullsEye.cpp \
    HayStack.cpp \
    Access.cpp \
    isl-helpers.cpp \
    Timer.cpp

libbullseye_la_CPPFLAGS = $(AM_CPPFLAGS) $(BOOST_CPPFLAGS) 
libbullseye_la_LIBADD = $(LIB_ISL) $(LIB_BARVINOK) $(LIB_PET) 

# build the bullseye program
bin_PROGRAMS = bullseye
bullseye_SOURCES = main.cpp 

bullseye_CPPFLAGS = $(AM_CPPFLAGS) $(BOOST_CPPFLAGS)
bullseye_LDFLAGS = -static $(BOOST_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIB) 
bullseye_LDADD = libbullseye.la $(LIB_ISL) $(LIB_BARVINOK) $(LIB_PET) $(CPLEX_HOME)cplex/lib/x86-64_linux/static_pic/libcplex.a $(CPLEX_HOME)cplex/lib/x86-64_linux/static_pic/libilocplex.a $(CPLEX_HOME)concert/lib/x86-64_linux/static_pic/libconcert.a $(CPLEX_HOME)opl/lib/x86-64_linux/static_pic/libilocplex.a $(CPLEX_HOME)/opl/lib/x86-64_linux/static_pic/libconcert.a 

