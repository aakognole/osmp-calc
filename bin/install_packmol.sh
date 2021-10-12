#!/bin/bash

if [ ! -e packmol/packmol ]; then
    wget https://github.com/m3g/packmol/archive/refs/heads/master.zip
    unzip -q master.zip
    rm -rf master.zip packmol
    mv packmol-master packmol
    cd packmol
    make
    cd ..
else
    exit
fi
