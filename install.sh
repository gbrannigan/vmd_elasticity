#!/bin/bash

PLUGINDIR=${HOME}/local/vmd-1.9.4/vmd/plugins/LINUXAMD64/tcl
VERSION=1.0

mkdir -p ${PLUGINDIR}/bending_modulus${VERSION}
rsync -av cmake-build-debug/libvmd_bending_modulus.so pkgIndex.tcl ${PLUGINDIR}/bending_modulus${VERSION}/. 
