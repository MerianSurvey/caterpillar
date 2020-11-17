This is the official release of the standard system response (instrument + atmosphere) for DES DR1 release (also known as Y3A2 Gold).
The system response is defined in steps of 5 Angstrom, for DES-grizY filter bandpasses, from 3800 Angstrom to 11000 Angstrom. The system responses (g/r/i/z/Y columns) includes the total throughput of the system (instrument+atmosphere). If you want the system throughput *without* the atmosphere, divide g/r/i/z/Y column by the atm column.

Instrument throughput is defined as focal plane average from DECal scan for in-band wavelengths, and is defined as zero for out-of-band wavelengths. The in-band wavelengths are defined as (in Angstrom):
g: 3800 - 5650
r: 5400 - 7350
i: 6750 - 8700
z: 8250 - 10150
Y: 9300 - 10700

Atmospheric throughput is computed using ModTran 4 with the following atmospheric parameters: 
Barometric Pressure 778 mb
Precipitable Water Vapor 3 mm
Ozone 263 Dobson
Aerosol Optical Depth 0.03
Aerosol Optical Index 1
Airmass 1.2
(Atmospheric throughput alone is available in the atm column)
