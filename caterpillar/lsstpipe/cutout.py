"""Generate postamps for Merian galaxies using HSC or DECam images.

Based on Johnny Greco's lsstutils package: https://github.com/johnnygreco/lsstutils
But upgraded to the Gen3 middleware

- 2022-02-17
"""

import os
import shutil
import argparse
import warnings

from datetime import date

import numpy as np

import astropy.units as u
from astropy.table import Table, QTable

from joblib import Parallel, delayed
from spherical_geometry.polygon import SphericalPolygon

try:
    import lsst.log
    Log = lsst.log.Log()
    Log.setLevel(lsst.log.ERROR)

    import lsst.daf.butler as dafButler

    import lsst.daf.base
    import lsst.geom as geom
    import lsst.afw.image as afwImage
    import lsst.afw.geom as afwGeom
    import lsst.afw.coord as afwCoord

    import lsst.afw.display
    import lsst.afw.display.rgb as afwRgb
except ImportError:
    warnings.warn("lsstPipe is not installed. Please install it first.")

MERIAN_REPO = '/projects/MERIAN/repo'
S18A_WIDE_ROOT = '/tigress/HSC/DR/s18a_wide'
PIXEL_SCALE = 0.168 # arcsec / pixel


def _prepare_dataset(root, collections):
    """
    Get the butler for the given dataset and return the skyMap object.
    """
    butler = dafButler.Butler(root, collections=collections)
    return butler, butler.get('skyMap')

def _get_ra_dec_name(id_arr, ra_arr, dec_arr):
    """Get the object name based on ID and (RA, Dec)."""
    return [
        "{:s}_{:s}_{:s}_{:s}".format(
            str(i), "{:8.4f}".format(ra).strip(), "{:8.4f}".format(dec).strip()
            ) for (i, ra, dec) in zip(id_arr, ra_arr, dec_arr)]

def _get_file_prefix(name_arr, band, prefix):
    """Get the prefix of the output files based on the ID."""
    if prefix is None:
        return ["{:s}_{:s}".format(str(name), band) for name in name_arr]
    else:
        return ["{:s}_{:s}_{:s}".format(prefix, str(name), band) for name in name_arr]

def _get_output_dir(output_dir, chunk_arr, name_arr):
    """Get the directory for the output cutout data."""
    # Check the output directory
    if not os.path.isdir(output_dir):
        raise ValueError("Output directory '{:s}' does not exist".format(output_dir))
    
    return [os.path.join(output_dir, str(chunk), str(name)) for (chunk, name) in zip(chunk_arr, name_arr)]
    
def _get_int_chunk(data, n_chunk):
    """Assign integer chunk ID to the data."""
    if n_chunk > len(data):
        raise ValueError("Too many chunks...")
    if n_chunk <= 0:
        raise ValueError("Chunk number has to be larger than 0...")
    
    chunk_arr = np.ones(len(data), dtype=int)
    if n_chunk == 1:
        return chunk_arr
    
    chunk_size = np.ceil(len(data) / n_chunk).astype(int)

    start, end = 0, chunk_size
    for i in np.arange(n_chunk):
        chunk_arr[start: end] = i + 1
        start, end = end, end + chunk_size 
        end = len(data) if end > len(data) else end
    
    return chunk_arr

def _prepare_input_cat(input_cat, half_size, unit, ra_col, dec_col, band, id_col, chunk, 
                       prefix, output_dir, save=True):
    """
    Prepare the input sample for the given dataset.
    
    The cutouts are organized into:
        [output_dir]/[chunk_id]/[galaxy_id]/[file_name].fits
    And the file name prefix is: 
        ([prefix]_[galaxy_id]_[band]
    """
    # Load the input catalog
    if isinstance(input_cat, str):
        input_cat = Table.read(input_cat)

    # Get an array for half size
    if isinstance(half_size, str):
        if half_size.strip() not in input_cat.colnames:
            raise ValueError("Wrong half size column name. [{:s}]".format(half_size))
        half_size_arr = input_cat[half_size]
    else:
        # Using the same size for all objects
        half_size_arr = np.full(len(input_cat), float(half_size))
    
    if np.any(half_size_arr < 0):
        raise ValueError("Negative size value.")
    
    # Add size unit if necessary
    if unit != 'pixel' and half_size_arr.unit is None:
        # Check the half size unit
        if unit.strip() not in ['arcsec', 'arcmin', 'degree', 'pixel']:
            raise ValueError("Wrong size unit. [arcsec, arcmin, degree, pixel]")
        half_size_arr = [s * u.Unit(unit) for s in half_size_arr]

    # Get the RA and DEC arrays
    if ra_col not in input_cat.colnames:
        raise ValueError("Wrong R.A. column name. [{:s}]".format(ra_col))
    if dec_col not in input_cat.colnames:
        raise ValueError("Wrong Dec column name. [{:s}]".format(dec_col))
    ra_arr, dec_arr = input_cat[ra_col], input_cat[dec_col]
    
    # Get the output directory and file name 
    
    # Get the object id or name
    if id_col is None:
        name_arr = _get_ra_dec_name(np.arange(len(ra_arr)) + 1, ra_arr, dec_arr)
    else:
        if id_col not in input_cat.colnames:
            raise ValueError("Wrong ID column name. [{:s}]".format(id_col))
        name_arr = input_cat[id_col]
    
    # Get the output file prefix 
    prefix_arr = _get_file_prefix(name_arr, band, prefix)
    
    # Get the directory of the output file
    if chunk is not None:
        if isinstance(chunk, str):
            if chunk not in input_cat.colnames:
                raise ValueError("Wrong Chunk column name. [{:s}]".format(chunk))
            chunk_arr = input_cat[chunk]
        else:
            chunk_arr = _get_int_chunk(input_cat, int(chunk))
    else:
        chunk_arr = None
        
    # Get the output file directory
    dir_arr = _get_output_dir(output_dir, chunk_arr, name_arr)

    sample = QTable(
        [name_arr, prefix_arr, dir_arr, chunk_arr, list(ra_arr), list(dec_arr), half_size_arr],
        names=('name', 'prefix', 'dir', 'chunk', 'ra', 'dec', 'half_size')
    )

    if save:
        today = date.today()
        prefix = 'postamps' if prefix is None else prefix
        sample.write(
            os.path.join(output_dir, "{:s}-{:4d}-{:02d}-{:02d}.fits".format(
                prefix, today.year, today.month, today.day)), overwrite=True)
    
    return sample

def _afw_coords(coord_list):
    """
    Convert list of ra and dec to lsst.afw.coord.IcrsCoord.

    Parameters
    ----------
    coord_list : list of tuples or tuple
        ra and dec in degrees.

    Returns
    -------
    afw_coords : list of lsst.afw.coord.IcrsCoord
    """
    if type(coord_list[0]) in (float, int, np.float64):
        ra, dec = coord_list
        afw_coords = geom.SpherePoint(ra * geom.degrees, dec * geom.degrees)
    else:
        afw_coords = [
            geom.SpherePoint(ra * geom.degrees, dec * geom.degrees) for ra, dec in coord_list]

    return afw_coords

def sky_cone(ra_c, dec_c, theta, steps=50, include_center=True):
    """
    Get ra and dec coordinates of a cone on the sky.

    Parameters
    ----------
    ra_c, dec_c: float
        Center of cone in degrees.
    theta: astropy Quantity, float, or int
        Angular radius of cone. Must be in arcsec
        if not a Quantity object.
    steps: int, optional
        Number of steps in the cone.
    include_center: bool, optional
        If True, include center point in cone.

    Returns
    -------
    ra, dec: ndarry
        Coordinates of cone.
    """
    if isinstance(theta, float) or isinstance(theta, int):
        theta = theta * u.Unit('arcsec')

    cone = SphericalPolygon.from_cone(
        ra_c, dec_c, theta.to('deg').value, steps=steps)
    ra, dec = list(cone.to_lonlat())[0]
    ra = np.mod(ra - 360., 360.0)
    if include_center:
        ra = np.concatenate([ra, [ra_c]])
        dec = np.concatenate([dec, [dec_c]])
    return ra, dec

def get_tract_patch_list(coord_list, skymap):
    """
    Find the tracts and patches that overlap with the
    coordinates in coord_list. Pass the four corners of
    a rectangle to get all tracts and patches that overlap
    with this region.

    Parameters
    ----------
    coord_list : list (tuples or lsst.afw.coord.IcrsCoord)
        ra and dec of region
    skymap : lsst.skymap.ringsSkyMap.RingsSkyMap, optional
        The lsst/hsc skymap.

    Returns
    -------
    region_ids : structured ndarray
        Tracts and patches that overlap coord_list.
    tract_patch_dict : dict
        Dictionary of dictionaries, which takes a tract
        and patch and returns a patch info object.
    """
    if isinstance(coord_list[0], float) or isinstance(coord_list[0], int):
        coord_list = [_afw_coords(coord_list)]
    elif not isinstance(coord_list[0], geom.SpherePoint):
        coord_list = _afw_coords(coord_list)

    tract_patch_list = skymap.findTractPatchList(coord_list)

    ids = []
    for tract_info, patch_info_list in tract_patch_list:
        for patch_info in patch_info_list:
            ids.append((tract_info.getId(), patch_info.getSequentialIndex()))

    return np.array(ids, dtype=[('tract', int), ('patch', int)])


def _get_patches(butler, skymap, coord_list, band, data_type='deepCoadd'):
    """
    Retrieve the data products for all the patches that overlap with the coordinate.
    """
    # Retrieve the Tracts and Patches that cover the cutout region
    patches = get_tract_patch_list(coord_list, skymap)

    # Collect the images
    images = []
    for t, p in patches:
        data_id = {'tract' : t, 'patch' : p, 'band' : band.upper()}
        if butler.datasetExists(data_type, data_id):
            img = butler.get(data_type, data_id, immediate=True)
            images.append(img)

    if len(images) == 0:
        return None
    return images

def _get_single_cutout(img, coord, half_size_pix):
    """Cutout from a single patch image.
    
    half_size_pix needs to be in pixels.
    """
    # Get the WCS and the pixel coordinate of the central pixel
    wcs = img.getWcs()
    pix = wcs.skyToPixel(coord)
    pix = geom.Point2I(pix)

    # Define a bounding box for the cutout region
    bbox = geom.Box2I(pix, pix)
    bbox.grow(half_size_pix)

    # Original pixel coordinate of the bounding box
    x0, y0 = bbox.getBegin()

    # Clip the cutout region from the original image
    bbox.clip(img.getBBox(afwImage.PARENT))

    # Make an afwImage object
    cut = img.Factory(img, bbox, afwImage.PARENT)

    return cut, x0, y0

def _build_cutout_wcs(coord, cutouts, index, origins):
    """Build new WCS header for the cutout."""
    # Get the WCS information from the largest cutout
    largest_cutout = cutouts[index]
    subwcs = largest_cutout.getWcs()

    # Information for the WCS header
    crpix_1, crpix_2 = subwcs.skyToPixel(coord)
    crpix_1 -= origins[index][0]
    crpix_2 -= origins[index][1]
    cdmat = subwcs.getCdMatrix()

    wcs_header = lsst.daf.base.PropertyList()
    wcs_header.add('CRVAL1', coord.getRa().asDegrees())
    wcs_header.add('CRVAL2', coord.getDec().asDegrees())
    wcs_header.add('CRPIX1', crpix_1 + 1)
    wcs_header.add('CRPIX2', crpix_2 + 1)
    wcs_header.add('CTYPE1', 'RA---TAN')
    wcs_header.add('CTYPE2', 'DEC--TAN')
    wcs_header.add('CD1_1', cdmat[0, 0])
    wcs_header.add('CD2_1', cdmat[1, 0])
    wcs_header.add('CD1_2', cdmat[0, 1])
    wcs_header.add('CD2_2', cdmat[1, 1])
    wcs_header.add('RADESYS', 'ICRS')

    return afwGeom.makeSkyWcs(wcs_header)

def _get_psf(exp, coord):
    """Get the coadd PSF image.

    Parameters
    ----------
    exp: lsst.afw.image.exposure.exposure.ExposureF
        Exposure
    coord: lsst.geom.SpherePoint
        Coordinate for extracting PSF

    Returns
    -------
    psf_img: lsst.afw.image.image.image.ImageD
        2-D PSF image
    """
    wcs = exp.getWcs()
    if not isinstance(coord, geom.SpherePoint):
        coord = _afw_coords(coord)
    coord = wcs.skyToPixel(coord)
    psf = exp.getPsf()

    try:
        psf_img = psf.computeKernelImage(coord)
        return psf_img
    except Exception:
        print('**** Cannot compute PSF Image *****')
        return None

def generate_cutout(butler, skymap, ra, dec, band='N708', data_type='deepCoadd',
                    half_size=10.0 * u.arcsec, psf=True, verbose=False):
    """Generate a single cutout image.
    """
    if not isinstance(half_size, u.Quantity):
        # Assume that this is in pixel
        half_size_pix = int(half_size)
    else:
        half_size_pix = int(half_size.to('arcsec').value / PIXEL_SCALE)

    # Width and height of the post-stamps
    stamp_shape = (half_size_pix * 2 + 1, half_size_pix * 2 + 1)

    # Coordinate of the image center
    coord = geom.SpherePoint(ra * geom.degrees, dec * geom.degrees)

    # Make a list of (RA, Dec) that covers the cutout region
    radec_list = np.array(
        sky_cone(ra, dec, half_size_pix * PIXEL_SCALE * u.Unit('arcsec'), steps=50)).T

    # Retrieve the Patches that cover the cutout region
    img_patches = _get_patches(butler, skymap, radec_list, band, data_type=data_type)

    if img_patches is None:
        if verbose:
            print('***** No data at {:.5f} {:.5f} *****'.format(ra, dec))
        return None

    # Making the stacked cutout
    cutouts = []
    idx, bbox_sizes, bbox_origins = [], [], []

    for img_p in img_patches:
        # Generate cutout
        cut, x0, y0 = _get_single_cutout(img_p, coord, half_size_pix)
        cutouts.append(cut)
        # Original lower corner pixel coordinate
        bbox_origins.append([x0, y0])
        # New lower corner pixel coordinate
        xnew, ynew = cut.getBBox().getBeginX() - x0, cut.getBBox().getBeginY() - y0
        idx.append([xnew, xnew + cut.getBBox().getWidth(),
                    ynew, ynew + cut.getBBox().getHeight()])
        # Area of the cutout region on this patch in unit of pixels
        # Will reverse rank all the overlapped images by this
        bbox_sizes.append(cut.getBBox().getWidth() * cut.getBBox().getHeight())

    # Stitch cutouts together with the largest bboxes inserted last
    stamp_bbox = geom.BoxI(geom.Point2I(0,0), geom.Extent2I(*stamp_shape))
    stamp = afwImage.MaskedImageF(stamp_bbox)
    bbox_sorted_ind = np.argsort(bbox_sizes)

    for i in bbox_sorted_ind:
        masked_img = cutouts[i].getMaskedImage()
        stamp[idx[i][0]: idx[i][1], idx[i][2]: idx[i][3]] = masked_img

    # Build the new WCS of the cutout
    stamp_wcs = _build_cutout_wcs(coord, cutouts, bbox_sorted_ind[-1], bbox_origins)

    cutout = afwImage.ExposureF(stamp, stamp_wcs)

    # The final product of the cutout
    if psf:
        psf = _get_psf(cutouts[bbox_sorted_ind[-1]], coord)
        return cutout, psf
    return cutout


def catalog_cutout(input_cat, root, collections, band, half_size=10, unit='arcsec', psf=True, data_type='deepCoadd_calexp',
                   output_dir='./', prefix=None, ra_col='ra', dec_col='dec', id_col='id', chunk=None):
    """
    Generate cutout/postamps for objects in a catalog using coadd images reduced by `lsstPipe`.

    Parameters
    ----------
    input_cat : `astropy.table.Table` or str
        Input catalog. Should contain at least the (RA, Dec) coordinates of the objects. 
    root : str
        Location of the data repo. 
    collections : str
        The collections sub-directory to the dataset.
    band : str
        Name of the filter.
    half_size : float or str
        Half size of the cutout image. 
        Could be the column name of the half size in the input catalog, or the half size value for the whole catalog.
        Used in combination with the unit.
    unit : str
        Unit of the half size.
    psf : bool
        Whether to generate PSF images.
    data_type : str
        The type of data to be used for the cutout.
    output_dir : str
        The output directory.
    prefix : str
        Prefix of the output file.
    ra_col : str
        Name of the RA column in the input catalog.
    dec_col : str
        Name of the Dec column in the input catalog.
    id_col : str
        Name of the ID column in the input catalog.
    chunk : int
        Number of chunks to split the catalog into.
    """

    # Check the required data type
    if data_type.strip() not in ['deepCoadd', 'deepCoadd_calexp', 'deepCoadd_background', 'deepCoadd_calexp_background']:
        raise ValueError("Wrong coadd type. [deepCoadd, deepCoadd_calexp, deepCoadd_background, deepCoadd_calexp_background]")

    # Get the butler for the given dataset and return the skyMap object.
    butler, skyMap = _prepare_dataset(root, collections)

    # Prepare the input catalog
    sample = _prepare_input_cat(
        input_cat, half_size, unit, ra_col, dec_col, band, id_col, chunk, prefix, output_dir, save=False)


def test_cutout():
    """
    Test the cutout function.
    """
    input_cat = '/home/sh19/work/cosmos-2022-02-17.fits'

    root = '/projects/MERIAN/repo'
    collections = 'DECam/runs/merian/w_2022_02/t9813_deep_N708'

    band = 'N708'
    data_type='deepCoadd_calexp'

    butler, skymap = _prepare_dataset(root, collections)

    # Prepare the input catalog
    sample = _prepare_input_cat(
        input_cat, 20.0, 'arcsec', 'ra', 'dec', 'N708', 'id', 20, 'cosmos', './', save=False)

    ra, dec = sample[0]['ra'], sample[0]['dec']
    print("RA = {:.5f}, Dec = {:.5f}".format(ra, dec))

    # Coordinate of the image center
    coord = geom.SpherePoint(ra * geom.degrees, dec * geom.degrees)

    # Make a list of (RA, Dec) that covers the cutout region
    radec_list = np.array(
        sky_cone(ra, dec, 100 * PIXEL_SCALE * u.Unit('arcsec'), steps=50)).T

    # Retrieve the Patches that cover the cutout region
    img_patches = _get_patches(butler, skymap, radec_list, band, data_type=data_type)

    return img_patches
