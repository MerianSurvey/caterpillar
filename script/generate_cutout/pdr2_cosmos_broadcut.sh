#!/bin/sh

# Setup hscPipe enviroment
module load rh/devtoolset/6
. /tigress/HSC/LSST/stack3_tiger/loadLSST.bash

setup lsst_apps
setup obs_subaru

export MERIAN="/tigress/MERIAN"

python3 ../merian_tiger_cutout.py \
    $MERIAN"/catalog/cosmos/pdr2_cosmos_wide_broadcut_bsm.fits" \
    --test cosmos_broad --chunk chunk_id \
    --ra ra --dec dec --name object_id \
    --prefix cosmos --size half_size --unit arcsec \
    --njobs 4 --psf
