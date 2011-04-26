AM_CFLAGS = $(DEBUG_FLAGS)
AM_LDFLAGS = -no-undefined -version-info 1:0:1

lib_LTLIBRARIES = libipl.la
libipl_la_LIBADD = $(LIBLEPT_LIBS) $(LIBFTTW3_LIBS) $(LIBFUZZYLITE_LIBS)

libipl_la_SOURCES = corner.c cornerlow.c dpix.c fourier.c fourierstub.c	    \
 frequency1.c frequency2.c frequency2stub.c motion.c motionlow.c noise.c    \
 pixafunc.c random.c spatial.c spatiallow.c spectra.c spectrastub.c         \
 window.c

pkginclude_HEADERS = ipl.h

LDADD = libipl.la
