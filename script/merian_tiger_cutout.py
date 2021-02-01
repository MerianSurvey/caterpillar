"""Generate HSC cutout for Merian survey on

Based on Johnny Greco's lsstutils package: https://github.com/johnnygreco/lsstutils
"""

import os
import shutil
import argparse

from datetime import date

import numpy as np

import astropy.units as u
from astropy.table import Table, QTable

from joblib import Parallel, delayed
from spherical_geometry.polygon import SphericalPolygon

import lsst.log
Log = lsst.log.Log()
Log.setLevel(lsst.log.ERROR)

import lsst.daf.base
import lsst.geom as geom
import lsst.daf.persistence as dafPersist
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.coord as afwCoord

import lsst.afw.display
import lsst.afw.display.rgb as afwRgb


MERIAN_ROOT = '/tigress/MERIAN/cutout'
DATA_ROOT = '/tigress/HSC/DR/s18a_wide'
PIXEL_SCALE = 0.168 # arcsec / pixel


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


def get_psf(exp, coord):
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
        coord = make_afw_coords(coord)
    coord = wcs.skyToPixel(coord)
    psf = exp.getPsf()

    try:
        psf_img = psf.computeKernelImage(coord)
        return psf_img
    except Exception:
        print('**** Cannot compute PSF Image *****')
        return None


def make_afw_coords(coord_list):
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
        afw_coords = geom.SpherePoint(
            geom.Angle(ra, geom.degrees),
            geom.Angle(dec, geom.degrees))
    else:
        afw_coords = [
            geom.SpherePoint(
                geom.Angle(ra, geom.degrees),
                geom.Angle(dec, geom.degrees)) for ra, dec in coord_list]

    return afw_coords


def tracts_n_patches(coord_list, skymap=None, data_dir=DATA_ROOT):
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
        The lsst/hsc skymap. If None, it will be created.
    data_dir : string, optional
        Rerun directory. Will use name in .superbutler
        by default.

    Returns
    -------
    region_ids : structured ndarray
        Tracts and patches that overlap coord_list.
    tract_patch_dict : dict
        Dictionary of dictionaries, which takes a tract
        and patch and returns a patch info object.
    """
    if isinstance(coord_list[0], float) or isinstance(coord_list[0], int):
        coord_list = [make_afw_coords(coord_list)]
    elif not isinstance(coord_list[0], geom.SpherePoint):
        coord_list = make_afw_coords(coord_list)

    if skymap is None:
        butler = lsst.daf.persistence.Butler(data_dir)
        skymap = butler.get('deepCoadd_skyMap', immediate=True)

    tract_patch_list = skymap.findTractPatchList(coord_list)

    ids = []
    tract_patch_dict = {}
    for tract_info, patch_info_list in tract_patch_list:
        patch_info_dict = {}
        for patch_info in patch_info_list:
            patch_index = patch_info.getIndex()
            patch_id = str(patch_index[0]) + ',' + str(patch_index[1])
            ids.append((tract_info.getId(), patch_id))
            patch_info_dict.update({patch_id : patch_info})
        tract_patch_dict.update({tract_info.getId() : patch_info_dict})

    region_ids = np.array(ids, dtype=[('tract', int), ('patch', 'S4')])

    return region_ids, tract_patch_dict

def make_single_cutout(img, coord, radius):
    """Cutout from a single patch image."""
    # Get the WCS and the pixel coordinate of the central pixel
    wcs = img.getWcs()
    pix = wcs.skyToPixel(coord)
    pix = geom.Point2I(pix)

    # Define a bounding box for the cutout region
    bbox = geom.Box2I(pix, pix)
    bbox.grow(radius)

    # Original pixel coordinate of the bounding box
    x0, y0 = bbox.getBegin()

    # Clip the cutout region from the original image
    bbox.clip(img.getBBox(afwImage.PARENT))

    # Make an afwImage object
    cut = img.Factory(img, bbox, afwImage.PARENT)

    return cut, x0, y0

def build_cutout_wcs(coord, cutouts, index, origins):
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

def generate_cutout(butler, skymap, ra, dec, band='i', label='deepCoadd_skyMap',
                    radius=10.0 * u.arcsec, psf=True, verbose=False):
    """Generate a single cutout image.
    """
    if not isinstance(radius, u.Quantity):
        # Assume that this is in pixel
        radius = int(radius)
    else:
        radius = int(radius.to('arcsec').value / PIXEL_SCALE)

    # Width and height of the post-stamps
    stamp_shape = (radius * 2 + 1, radius * 2 + 1)

    # Coordinate of the image center
    coord = geom.SpherePoint(
        ra * geom.degrees, dec * geom.degrees)

    # Make a list of (RA, Dec) that covers the cutout region
    radec_list = np.array(sky_cone(ra, dec, radius * PIXEL_SCALE, steps=50)).T

    # Retrieve the Tracts and Patches that cover the cutout region
    patches, _ = tracts_n_patches(radec_list, skymap)

    # Collect the images
    images = []
    for t, p in patches:
        data_id = {'tract' : t, 'patch' : p.decode(),
                   'filter' : 'HSC-' + band.upper()}

        if butler.datasetExists(label, data_id):
            img = butler.get(label, data_id, immediate=True)
            images.append(img)

    if len(images) == 0:
        if verbose:
            print('***** No data at {:.5f} {:.5f} *****'.format(ra, dec))
        return None

    cutouts = []
    idx, bbox_sizes, bbox_origins = [], [], []

    for img_patch in images:
        # Generate cutout
        cut, x0, y0 = make_single_cutout(img_patch, coord, radius)
        cutouts.append(cut)

        # Original lower corner pixel coordinate
        bbox_origins.append([x0, y0])

        # New lower corner pixel coordinate
        xnew, ynew = cut.getBBox().getBeginX() - x0, cut.getBBox().getBeginY() - y0
        idx.append([xnew, xnew + cut.getBBox().getWidth(),
                    ynew, ynew + cut.getBBox().getHeight()])

        # Pixel size of the cutout region
        bbox_sizes.append(cut.getBBox().getWidth() * cut.getBBox().getHeight())

    # Stitch cutouts together with the largest bboxes inserted last
    stamp_bbox = geom.BoxI(geom.Point2I(0,0), geom.Extent2I(*stamp_shape))
    stamp = afwImage.MaskedImageF(stamp_bbox)
    bbox_sorted_ind = np.argsort(bbox_sizes)

    for i in bbox_sorted_ind:
        masked_img = cutouts[i].getMaskedImage()
        stamp[idx[i][0]: idx[i][1], idx[i][2]: idx[i][3]] = masked_img

    # Build the new WCS of the cutout
    stamp_wcs = build_cutout_wcs(coord, cutouts, bbox_sorted_ind[-1], bbox_origins)

    # The final product of the cutout
    if psf:
        return (afwImage.ExposureF(stamp, stamp_wcs),
                get_psf(cutouts[bbox_sorted_ind[-1]], coord))
    return afwImage.ExposureF(stamp, stamp_wcs)

def cutout_one(butler, skymap, obj, band, label, psf):
    """Generate cutout for a single object."""
    prefix, ra, dec, radius = obj['prefix'], obj['ra'], obj['dec'], obj['radius']

    cutout = generate_cutout(
        butler, skymap, ra, dec, band=band, label=label, radius=radius, psf=psf)

    # Make a new folder is necessary
    if not os.path.isdir(os.path.split(prefix)[0]):
        os.makedirs(os.path.split(prefix)[0], exist_ok=True)

    if psf:
        img, psf = cutout
        img.writeFits("{:s}.fits".format(prefix))
        psf.writeFits("{:s}_psf.fits".format(prefix))
    else:
        cutout.writeFits("{:s}.fits".format(prefix))

def is_number(string):
    """Check if a string can be converted into a float number."""
    try:
        float(string)
        return True
    except ValueError:
        return False

def ra_dec_prefix(ra_arr, dec_arr, band, prefix):
    """Get the output file name based on (RA, Dec)."""
    return [
        "{:s}_{:s}_{:s}_{:s}".format(
            prefix, "{:8.4f}".format(ra).strip(), "{:8.4f}".format(dec).strip(),
            band.lower().strip()
            ) for (ra, dec) in zip(ra_arr, dec_arr)]

def prepare_catalog_merian(cat, size, band, ra='ra', dec='dec', name=None, unit='arcsec',
                           prefix=None, merian_root=MERIAN_ROOT, rerun=None, chunk=None):
    """Get the input catalog ready for generating cutouts for Merian survey.

    The format for data storage of Merian cutous on Tiger is:
        MERIAN/cutout/[RERUN]/[CHUNK]/[GALAXY_ID]/[DATA_SET]

    Here:
        - `DATA_SET` is "hsc"
    """
    # RA, Dec: convert to standard column names
    ra_arr, dec_arr = cat[ra], cat[dec]

    # Folder for the dataset
    if rerun is None:
        rerun = 'test'
    merian_output = os.path.join(merian_root, rerun)

    if not os.path.isdir(merian_output):
        os.mkdir(merian_output, exist_ok=True)

    if chunk is None:
        chunk_arr = np.full(len(ra_arr), 1)
    else:
        chunk_arr = cat[chunk]

    # Name of the output file
    if prefix is None:
        prefix = DATA_ROOT.split('/')[-1]
    prefix = prefix.lower()

    if name is None:
        id_arr = np.arange(len(ra_arr)) + 1
        name_arr = ra_dec_prefix(ra_arr, dec_arr, band, prefix)
    else:
        id_arr = cat[name]
        name_arr = [
            "{:s}_{:s}_{:s}".format(prefix, str(name).strip(), band.lower()) for name in cat[name]]

    output_arr = [
        os.path.join(merian_output, str(chunk).strip(), str(ids).strip(), name.strip())
        for (chunk, ids, name) in zip(chunk_arr, id_arr, name_arr)
    ]

    # Radius of the cutout
    if is_number(size):
        # Using the same size for all objects
        size_arr = np.full(len(cat), float(size))
    else:
        size_arr = cat[size]

    # Add size unit if necessary
    if unit != 'pixel':
        size_arr = [s * u.Unit(unit) for s in size_arr]

    sample = QTable(
        [name_arr, output_arr, list(ra_arr), list(dec_arr), size_arr],
        names=('name', 'prefix', 'ra', 'dec', 'radius')
    )

    today = date.today()
    sample.write(
        os.path.join(merian_output, "{:s}-{:4d}-{:02d}-{:02d}.fits".format(
            rerun, today.year, today.month, today.day)), overwrite=True)

    return sample


def batch_cutout_merian(incat, size=10, ra='ra', dec='dec', name=None, prefix=None,
                        unit='arcsec', label='deepCoadd_calexp', root=DATA_ROOT,
                        merian_root=MERIAN_ROOT, rerun=None, chunk=None,
                        njobs=1, psf=True):
    """Run the batch cutout script."""

    # Get the butler and skymap of the HSC data rerun
    butler = dafPersist.Butler(root)
    skymap = butler.get('deepCoadd_skyMap', immediate=True)

    if label.strip() not in ['deepCoadd', 'deepCoadd_calexp']:
        raise ValueError("Wrong coadd type. [deepCoadd, deepCoadd_calexp]")

    if unit.strip() not in ['arcsec', 'arcmin', 'degree', 'pixel']:
        raise ValueError("Wrong size unit. [arcsec, arcmin, degree, pixel]")

    # Read the input catalog
    cat = Table.read(incat)

    # Folder for the dataset
    if rerun is None:
        rerun = 'test'
    merian_output = os.path.join(merian_root, rerun)
    if not os.path.isdir(merian_output):
        os.makedirs(merian_output, exist_ok=True)

    for band in ['g', 'r', 'i', 'z', 'y']:
        print("# Will generate {:d} cutouts in {:s} band".format(len(cat), band))
        # Get the (RA, Dec)
        input_cat = prepare_catalog_merian(
            cat, size, band, ra=ra, dec=dec, name=name, unit=unit, prefix=prefix,
            merian_root=merian_root, rerun=rerun, chunk=chunk)

        if njobs <= 1:
            _ = [cutout_one(
                butler, skymap, obj, band, label, psf) for obj in input_cat]
        else:
            Parallel(n_jobs=njobs)(
                delayed(cutout_one)(
                    butler, skymap, obj, band, label, psf) for obj in input_cat)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input catalog")
    parser.add_argument(
        '-r', '--root', dest='root', help="Root directory of data repository",
        default=DATA_ROOT)
    parser.add_argument(
        '-l', '--label', dest='label', help="Label of the Coadd data product",
        default='deepCoadd_calexp')
    parser.add_argument(
        '-m', '--merian', dest='merian_root', help="Output directory",
        default=MERIAN_ROOT)
    parser.add_argument(
        '-t', '--test', dest='rerun', help="Name of the test or rerun",
        default=None)
    parser.add_argument(
        '-c', '--chunk', dest='chunk', help="Column for Tract/Brick or any other chunk ID",
        default=None)
    parser.add_argument(
        '-ra', '--ra', dest='ra', help="Column name for RA", default='ra')
    parser.add_argument(
        '-dec', '--dec', dest='dec', help="Column name for Dec", default='dec')
    parser.add_argument(
        '-n', '--name', dest='name', help="Column name for galaxy name or index",
        default=None)
    parser.add_argument(
        '-p', '--prefix', dest='prefix', help='Prefix of the output file',
        default=None)
    parser.add_argument(
        '-s', '--size', dest='size', help="Column name for the cutout size",
        default='radius')
    parser.add_argument(
        '-u', '--unit', dest='size_unit', help="Unit of the size column",
        default='arcsec')
    parser.add_argument(
        '-j', '--njobs', type=int, help='Number of jobs run at the same time',
        dest='njobs', default=1)
    parser.add_argument(
        '-psf', '--psf', action="store_true", dest='psf', default=False)

    args = parser.parse_args()

    batch_cutout_merian(
        args.input, size=args.size, band=args.band, ra=args.ra, dec=args.dec,
        name=args.name, prefix=args.prefix, unit=args.size_unit, merian_root=args.merian_root,
        chunk=args.chunk, rerun=args.rerun, label=args.label, root=args.root,
        njobs=args.njobs, psf=args.psf)
