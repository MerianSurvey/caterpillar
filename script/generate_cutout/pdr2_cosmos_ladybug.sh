#!/bin/sh

# Setup hscPipe enviroment
module load rh/devtoolset/6
. /tigress/HSC/LSST/stack3_tiger/loadLSST.bash

setup lsst_apps
setup obs_subaru

export MERIAN="/tigress/MERIAN"

python3 ../merian_tiger_cutout.py \
    $MERIAN"/catalog/cosmos/ladybug/ladybug_hsc_size.fits" \
    --test ladybug \
    --ra ra_cosmos --dec dec_cosmos --name index \
    --prefix ladybug --size half_size --unit arcsec \
    --njobs 1 --psf

python3 ../merian_tiger_cutout.py \
    $MERIAN"/catalog/cosmos/ladybug/ladybug_cl_hsc_size.fits" \
    --test ladybug_cl \
    --ra ra_cosmos --dec dec_cosmos --name index \
    --prefix ladybug_cl --size half_size --unit arcsec \
    --njobs 1 --psf

python3 ../merian_tiger_cutout.py \
    $MERIAN"/catalog/cosmos/ladybug/ladybug_ch_hsc_size.fits" \
    --test ladybug_ch \
    --ra ra_cosmos --dec dec_cosmos --name index \
    --prefix ladybug_ch --size half_size --unit arcsec \
    --njobs 1 --psf
