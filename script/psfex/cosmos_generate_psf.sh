#!/bin/sh
#Building an initial model of the PSF with with SExtractor and PSFEx:
sex /Users/yifei/Dropbox/merian/Merian_data/COSMOS/wide/c4d_210306_011257_osj_N708_deep.fits -c cosmos_deep_init.sex
psfex cosmos_deep_init.fits -c cosmos_deep.psfex
#Re-run the SExtractor with the PSF model
sex -PSF_NAME cosmos_deep_init.psf /Users/yifei/Dropbox/merian/Merian_data/COSMOS/wide/c4d_210306_011257_osj_N708_deep.fits -c cosmos_deep.sex
