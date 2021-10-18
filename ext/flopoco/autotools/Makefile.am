# Makefile for flopoco

SUBDIRS = src
DIST_SUBDIRS = src
ACLOCAL_AMFLAGS = -I m4

## Comment out to make new version of FloPoColib
## TODO: emulate CMakefiles.txt to link with Factories.cpp etc.
bin_PROGRAMS = flopoco longacc2fp fp2bin bin2fp fpadder_example




flopoco_SOURCES = src/main.cpp
flopoco_CPPFLAGS = $(CPPFLAGS) -I$(top_srcdir)/src/
flopoco_CXXFLAGS = $(CXXFLAGS)
flopoco_LDADD = $(LIBS) src/libflopoco.la




bin2fp_SOURCES = \
                 src/Tools/bin2fp.cpp \
                 src/utils.cpp
bin2fp_CPPFLAGS = $(CPPFLAGS) -I$(top_srcdir)/src/
bin2fp_CXXFLAGS = $(CXXFLAGS)
bin2fp_LDADD = $(LIBS) src/libflopoco.la




fp2bin_SOURCES = \
                 src/Tools/fp2bin.cpp \
                 src/utils.cpp
fp2bin_CPPFLAGS = $(CPPFLAGS) -I$(top_srcdir)/src/
fp2bin_CXXFLAGS = $(CXXFLAGS)
fp2bin_LDADD = $(LIBS) src/libflopoco.la



longacc2fp_SOURCES = \
                     src/Tools/longacc2fp.cpp \
                     src/utils.cpp
longacc2fp_CPPFLAGS = $(CPPFLAGS) -I$(top_srcdir)/src/
longacc2fp_CXXFLAGS = $(CXXFLAGS)
longacc2fp_LDADD = $(LIBS) src/libflopoco.la




fpadder_example_SOURCES = src/main_minimal.cpp
fpadder_example_CPPFLAGS = $(CPPFLAGS) -I$(top_srcdir)/src/
fpadder_example_CXXFLAGS = $(CXXFLAGS)
fpadder_example_LDADD = $(LIBS) src/libflopoco.la
