
VERSION=1.0
TCLINC=/usr/include/tcl
# Change this to point to your local VMD installation!
PLUGINDIR=${HOME}/local/vmd-1.9.3/plugins/LINUXAMD64/tcl

CPP=g++
#CPPFLAGS=-fpic -g -I${TCLINC}
CPPFLAGS=-fpic -O3 -I${TCLINC}
# LIBS=-Wl,-Bstatic -lfftw -Wl,-Bdynamic
LIBS=-lfftw3f

all: old_vmd_bending_modulus.so bending_modulus bending_modulus.tar.gz

bending_modulus: bending_modulus.cpp floatarray2d.h
	$(CPP) -g -O3 bending_modulus.cpp -o bending_modulus $(LIBS)

old_vmd_bending_modulus.so: old_vmd_bending_modulus.o
	$(CPP) -shared old_vmd_bending_modulus.o -o old_vmd_bending_modulus.so $(LIBS)

bending_modulus.tar.gz: old_vmd_bending_modulus.cpp Makefile
	tar czf bending_modulus${VERSION}.tar.gz old_vmd_bending_modulus.cpp Makefile pkgIndex.tcl

clean:
	rm *.o *.so *.tar.gz

install: all
	mkdir -p ${PLUGINDIR}/bending_modulus${VERSION} && \
		rsync -av old_vmd_bending_modulus.so pkgIndex.tcl ${PLUGINDIR}/bending_modulus${VERSION}/. 
