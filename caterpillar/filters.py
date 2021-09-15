#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Transmission curve related."""

__all__ = ['make_kfilter_file']


def make_kfilter_file(wave, curve, par_name='merian_filter.par'):
    """Save a transmission curve into a KFILTER file used by sedpy and kcorrect."""
    assert len(wave) == len(curve), "Something wrong with the curve!"
    
    with open(par_name, 'w') as par_file:
        # Print out the header
        par_file.write('typedef struct {\n')
        par_file.write('  double lambda;\n')
        par_file.write('  double pass;\n')
        par_file.write('} KFILTER; \n \n')
        
        for w, p in zip(wave, curve):
            par_file.write("  KFILTER    {:6.1f}    {:f}\n".format(w, p))
    
    return
