# Makefile for flopoco

SUBDIRS = src
DIST_SUBDIRS = src
ACLOCAL_AMFLAGS = -I m4


bin_PROGRAMS = flopoco longacc2fp fp2bin bin2fp ieee2bin bin2ieee fpadder_example




flopoco_SOURCES = \
                 src/main.cpp
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




bin2ieee_SOURCES = \
                 src/Tools/bin2ieee.cpp \
                 src/utils.cpp
bin2ieee_CPPFLAGS = $(CPPFLAGS) -I$(top_srcdir)/src/
bin2ieee_CXXFLAGS = $(CXXFLAGS)
bin2ieee_LDADD = $(LIBS) src/libflopoco.la




ieee2bin_SOURCES = \
                 src/Tools/ieee2bin.cpp \
                 src/utils.cpp
ieee2bin_CPPFLAGS = $(CPPFLAGS) -I$(top_srcdir)/src/
ieee2bin_CXXFLAGS = $(CXXFLAGS)
ieee2bin_LDADD = $(LIBS) src/libflopoco.la




longacc2fp_SOURCES = \
                     src/Tools/longacc2fp.cpp \
                     src/utils.cpp
longacc2fp_CPPFLAGS = $(CPPFLAGS) -I$(top_srcdir)/src/
longacc2fp_CXXFLAGS = $(CXXFLAGS)
longacc2fp_LDADD = $(LIBS) src/libflopoco.la




fpadder_example_SOURCES = src/main_minimal.cpp src/utils.cpp
fpadder_example_CPPFLAGS = $(CPPFLAGS) -I$(top_srcdir)/src/
fpadder_example_CXXFLAGS = $(CXXFLAGS)
fpadder_example_LDADD = $(LIBS) src/libflopoco.la
