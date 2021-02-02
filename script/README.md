# Python Scripts for Merian Data Reduction and Analysis

- Last Update: 2021-02-02

----

### Generating HSC cutouts on `Tiger`: `merian_tiger_cutout.py`

- This script will take an input FITS catalog of target galaxies and generate five-band HSC cutout images (flux, inverse variance, and mask) and PSF model images. It will also store the cutout data in a pre-designed folder structure:
    - Current structure: `MERIAN/cutout/[TEST_NAME]/[CHUNK_ID]/[GALAXY_ID]/[DATA_SET]`
    - `TEST_NAME`: separates different input catalog or different iteration of data reduction tests.
    - `CHUNK_ID`: separate large catalogs into small chunks. This could be HSC `Tract` ID or DECaLS `brick` or `sweep` ID.
    - `GALAXY_ID`: is the name or index of the galaxy in the catalog.
    - `DATA_SET`: separates imaging data from different sources. Here it should be `hsc`.
- The script will also save a copy of the input catalog with the minimum amount of useful information to the `TEST_NAME` folder. This includes the coordinates and the names of galaxies in the input catalog along with the size of the cutout.


- This requires local archive of HSC data (currently, only `S18A` is available on `Tiger`) and the LSST DM stack.
    - Current HSC data are available at `/tigress/HSC/DR/s18a_wide`
    - Current Merian data are stored at `/tigress/MERIAN/cutout`

#### Example Usage:

```bash
python3 merian_tiger_cutout.py ../data/cosmos_dwarf_test_10_2021Feb01.fits -t cosmos_test_1 -n id -s 20 -u arcsec --psf --ra ra --dec dec 
``` 

- Input catalog: `../data/cosmos_dwarf_test_10_2021Feb01.fits`
- `--ra`, `--dec`: Name of the columns for the coordinates of the galaxies. Here are `ra`, `dec`.
- `-t`: Test name is `cosmos_test_1`.
    - When not provided, the script will make a `test` folder for the catalog.
- `-n`: Column name for galaxy ID or name. Here is `id`
    - If not available, the script will name each object using its index in the input catalog.
- `-s`: Information about cutout size.
    - Here the value is a single number, which means this size will be applied to all galaxies.
    - If it is a string, like `size`, `radius`, then the script will look for this column in the catalog. This is useful for generating cutouts of different sizes for different galaxies.
- `-s`: Unit of the size. Here it is `arcsec`. 
    - Option includes: `arcsec`, `arcmin`, `degree`, `pixel`.
- `--psf`: Flag to turn on the PSF model generation. Without it, the script will only generate cutout images.


