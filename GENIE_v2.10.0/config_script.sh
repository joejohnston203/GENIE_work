#!/bin/bash

# Script to configure GENIE before building and installing

export GENIE="/home/joe/GENIE_work/GENIE_v2.10.0"
export ROOTSYS="/opt/root"

./configure \
    --prefix=/home/joe/GENIE_work \
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
    --with-pythia6-lib=/opt/v6_424/lib \
    --with-lhapdf-inc=/opt/lhapdf-5.9.1/include \
    --with-lhapdf-lib=/opt/lhapdf-5.9.1/lib \
    --with-libxml2-inc=/usr/include/libxml2 \
    --with-libxml2-lib=/usr/share/aclocal \
    --with-log4cpp-inc=/opt/log4cpp/include \
    --with-log4cpp-lib=/opt/log4cpp/lib
