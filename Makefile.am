if BUNDLED_ISL
    MAYBE_ISL = isl
endif
if BUNDLED_BARVINOK
    MAYBE_BARVINOK = barvinok
endif
if BUNDLED_PET
    MAYBE_PET = pet
endif

SUBDIRS = $(MAYBE_ISL) $(MAYBE_BARVINOK) $(MAYBE_PET) tests

ACLOCAL_AMFLAGS = -I m4

# -------------------------------------------------------------

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

libisl.la:
	$(MAKE) $(AM_MAKEFLAGS) -C isl libisl.la
libbarvinok.la:
	$(MAKE) $(AM_MAKEFLAGS) -C barvinok libbarvinok.la
libpet.la: 
	$(MAKE) $(AM_MAKEFLAGS) -C pet libpet.la

# -------------------------------------------------------------


MAINTAINERCLEANFILES = Makefile.in

AM_CPPFLAGS = @ISL_CFLAGS@ @BARVINOK_CFLAGS@ @PET_CFLAGS@ -I$(top_builddir)/
AM_LDFLAGS =  -lcplex -lconcert -lilocplex -lm -lpthread -ldl -lstdc++

LIB_ISL = $(ISL_LA) @ISL_LIBS@
LIB_BARVINOK = $(BARVINOK_LA) @BARVINOK_LIBS@
LIB_PET = $(PET_LA) @PET_LIBS@
LIB_CPLEX = $(CPLEX_LA) @CPLEX_LIBS@
LINK = $(CXXLINK)

lib_LTLIBRARIES = libbullseye.la
libbullseye_la_DEPENDENCIES = libpet.la libbarvinok.la libisl.la
libbullseye_la_HEADERS = \
    lib/Access.h \
    lib/BullsEye.h \
    lib/bullseyelib.h \
    lib/Definitions.h\
    lib/HayStack.h \
    lib/isl-helpers.h \
    lib/op.h \
    lib/Program.h \
    lib/Timer.h

libbullseye_la_SOURCES = \
    lib/Access.cpp \
    lib/Access.h \
    lib/BullsEye.cpp \
    lib/BullsEye.h \
    lib/bullseyelib.cpp \
    lib/bullseyelib.h \
    lib/Definitions.h\
    lib/HayStack.cpp \
    lib/HayStack.h \
    lib/isl-helpers.cpp \
    lib/isl-helpers.h \
    lib/op.h \
    lib/Program.cpp \
    lib/Program.h \
    lib/Timer.cpp \
    lib/Timer.h



libbullseye_la_CPPFLAGS = $(AM_CPPFLAGS) $(BOOST_CPPFLAGS) 
libbullseye_la_LDFLAGS = -static
libbullseye_la_LIBADD = $(LIB_ISL) \
                        $(LIB_BARVINOK) \
                        $(LIB_PET) \
                        $(CPLEX_HOME)cplex/lib/x86-64_linux/static_pic/libcplex.a \
                        $(CPLEX_HOME)cplex/lib/x86-64_linux/static_pic/libilocplex.a \
                        $(CPLEX_HOME)concert/lib/x86-64_linux/static_pic/libconcert.a \
                        $(CPLEX_HOME)opl/lib/x86-64_linux/static_pic/libilocplex.a \
                        $(CPLEX_HOME)/opl/lib/x86-64_linux/static_pic/libconcert.a \
                        $(AM_LDFLAGS)

# build the bullseye program
bin_PROGRAMS = bullseye
bullseye_DEPENDENCIES = libbullseye.la
bullseye_SOURCES = src/main.cpp libbullseye.la
bullseye_CPPFLAGS = $(AM_CPPFLAGS) $(BOOST_CPPFLAGS)
bullseye_LDFLAGS = -static $(BOOST_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIB) 
bullseye_LDADD = libbullseye.la $(LIB_ISL) $(LIB_BARVINOK) $(LIB_PET) $(CPLEX_HOME)cplex/lib/x86-64_linux/static_pic/libcplex.a $(CPLEX_HOME)cplex/lib/x86-64_linux/static_pic/libilocplex.a $(CPLEX_HOME)concert/lib/x86-64_linux/static_pic/libconcert.a $(CPLEX_HOME)opl/lib/x86-64_linux/static_pic/libilocplex.a $(CPLEX_HOME)/opl/lib/x86-64_linux/static_pic/libconcert.a 


libbullseye_ladir = $(includedir)/bullseye
include_HEADERS = $(libbullseye_la_HEADERS)














#BULLSEYE_LA = $(top_builddir)/lib/libbullseye.la
#AM_CPPFLAGS = @ISL_CFLAGS@ @BARVINOK_CFLAGS@ @PET_CFLAGS@  -I$(top_builddir)/
#AM_LDFLAGS =  -lcplex -lconcert -lilocplex -lm -lpthread -ldl -lstdc++ 
#LIB_ISL = $(ISL_LA) @ISL_LIBS@
#LIB_BARVINOK = $(BARVINOK_LA) @BARVINOK_LIBS@
#LIB_PET = $(PET_LA) @PET_LIBS@
#LIB_CPLEX = $(CPLEX_LA) @CPLEX_LIBS@
#LINK = $(CXXLINK)

# build the bullseye library
#lib_LIBRARIES = libbullseye.a
#libbullseye_la_SOURCES = \
#    lib/bullseyelib.h \
#    lib/bullseyelib.cpp


#ACLOCAL_AMFLAGS = -I m4
#
#bin_PROGRAMS = code
#lib_LIBRARIES = libcode.a
#lib_LTLIBRARIES = libcodes.la
#
#code_SOURCES = code.c libo.h libp.h
#code_LDFLAGS = -rpath $(libdir)
#code_LDADD = libcodes.la libcode.a
## or LDADD = libcodes.la libcode.a
#
#libcode_a_SOURCES = functionsp.c
#
#libcodes_la_SOURCES = functionso.c
#libcodes_la_LDFLAGS = -version-info 1:0:0


#bullseye: libbullseye.la
#	$(MAKE) $(MFLAGS) -C src bullseye

#info-bullseye:
#	$(MAKE) $(MFLAGS) -C src info-bullseye

#test_libbullseye: tests/LibbullseyeALinkTest.cpp
#    libtool --mode=link --tag=CC $(CC)  $
