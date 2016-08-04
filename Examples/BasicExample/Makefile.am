
AM_CPPFLAGS = -I${top_srcdir} -I${top_builddir}

if AUTO_EXAMPLES
bin_PROGRAMS = BasicExample
TESTS = BasicExample$(EXEEXT)
else
noinst_PROGRAMS = BasicExample
endif
SwapValuation_SOURCES = basicexample.cpp
SwapValuation_LDADD = ../../ql/libQuantLib.la ${BOOST_THREAD_LIB}

EXTRA_DIST = \
    CMakeLists.txt \
    BasicExample.vcxproj \
    BasicExample.vcxproj.filters \
    README.txt

.PHONY: examples check-examples

examples: BasicExample$(EXEEXT)

check-examples: examples
	./BasicExample$(EXEEXT)

dist-hook:
	mkdir -p $(distdir)/bin
	mkdir -p $(distdir)/build

