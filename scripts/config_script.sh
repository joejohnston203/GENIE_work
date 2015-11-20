#!/bin/bash

export GENIE="/home/joe/work/GENIE_v2.10.0"
export ROOTSYS="/usr/local/root"

./configure \
    --prefix=/usr \
    --disable-profiler \
    --enable-validation-tools \
    --disable-neut-cascade \
    --disable-cernlib  \
    --enable-lhapdf \
    --enable-flux-drivers \
    --enable-geom-drivers \
    --disable-doxygen \
    --enable-test \
    --enable-mueloss \
    --enable-dylibversion \
    --enable-event-server \
    --enable-t2k \
    --enable-numi \
    --enable-atmo \
    --enable-nucleon-decay \
    --enable-rwght \
    --enable-masterclass \
    --disable-debug \
    --with-optimiz-level=O2 \
    --with-pythia6-lib=/opt/pythia/v6_424/lib \
    --with-lhapdf-inc=/opt/LHAPDF-6.1.5/include \
    --with-lhapdf-lib=/opt/LHAPDF-6.1.5/src/.libs \
    --with-libxml2-inc=/usr/include/libxml2 \
    --with-libxml2-lib=/usr/share/aclocal \
    --with-log4cpp-inc=/opt/log4cpp/include \
    --with-log4cpp-lib=/opt/log4cpp/lib
