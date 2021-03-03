SELECT
	-- Basic information
	f1.object_id, f1.parent_id, f1.ra, f1.dec, f1.tract, f1.patch,
	-- Galactic extinction for correction
	f1.a_g, f1.a_r, f1.a_i, f1.a_z, f1.a_y,
	-- CModel photometry
	-- 1. Exponential component photometry (useful for dwarfs)
	f1.g_cmodel_exp_flux, f1.g_cmodel_exp_fluxsigma,
	f1.r_cmodel_exp_flux, f1.r_cmodel_exp_fluxsigma,
	f1.i_cmodel_exp_flux, f1.i_cmodel_exp_fluxsigma,
	f1.z_cmodel_exp_flux, f1.z_cmodel_exp_fluxsigma,
	f1.y_cmodel_exp_flux, f1.y_cmodel_exp_fluxsigma,
	-- 2. CModel photometry
	f1.g_cmodel_flux, f1.g_cmodel_fluxsigma,
	f1.r_cmodel_flux, f1.r_cmodel_fluxsigma,
	f1.i_cmodel_flux, f1.i_cmodel_fluxsigma,
	f1.z_cmodel_flux, f1.z_cmodel_fluxsigma,
	f1.y_cmodel_flux, f1.y_cmodel_fluxsigma,
	f1.g_cmodel_mag, f1.g_cmodel_magsigma,
	f1.r_cmodel_mag, f1.r_cmodel_magsigma,
	f1.i_cmodel_mag, f1.i_cmodel_magsigma,
	f1.z_cmodel_mag, f1.z_cmodel_magsigma,
	f1.y_cmodel_mag, f1.y_cmodel_magsigma,
	-- 3. Flags for CModel photometry
	-- If flag is set to True, it means the photometry is not reliable
	f1.g_cmodel_flag,
	f1.r_cmodel_flag,
	f1.i_cmodel_flag,
	f1.z_cmodel_flag,
	f1.y_cmodel_flag,
	-- PSF photometry
	f2.g_psfflux_flux as g_psf_flux, f2.g_psfflux_fluxsigma as g_psf_fluxsigma,
	f2.r_psfflux_flux as r_psf_flux, f2.r_psfflux_fluxsigma as r_psf_fluxsigma,
	f2.i_psfflux_flux as i_psf_flux, f2.i_psfflux_fluxsigma as i_psf_fluxsigma,
	f2.z_psfflux_flux as z_psf_flux, f2.z_psfflux_fluxsigma as z_psf_fluxsigma,
	f2.y_psfflux_flux as y_psf_flux, f2.y_psfflux_fluxsigma as y_psf_fluxsigma,
	-- PSF photometry flag (for quality cut later)
	f2.g_psfflux_flag as g_psf_flag,
	f2.r_psfflux_flag as r_psf_flag,
	f2.i_psfflux_flag as i_psf_flag,
	f2.z_psfflux_flag as z_psf_flag,
	f2.y_psfflux_flag as y_psf_flag,
	-- PSF-corrected aperture photometry
	-- 2_15: 1.1 arcsec seeing; 1.5 arcsec diamter aperture
	f4.g_convolvedflux_2_15_flux, f4.g_convolvedflux_2_15_fluxsigma,
	f4.r_convolvedflux_2_15_flux, f4.r_convolvedflux_2_15_fluxsigma,
	f4.i_convolvedflux_2_15_flux, f4.i_convolvedflux_2_15_fluxsigma,
	f4.z_convolvedflux_2_15_flux, f4.z_convolvedflux_2_15_fluxsigma,
	f4.y_convolvedflux_2_15_flux, f4.y_convolvedflux_2_15_fluxsigma,
	f4.g_convolvedflux_2_15_flag,
	f4.r_convolvedflux_2_15_flag,
	f4.i_convolvedflux_2_15_flag,
	f4.z_convolvedflux_2_15_flag,
	f4.y_convolvedflux_2_15_flag,
	-- 3_20: 1.3 arcsec seeing; 2.0 arcsec diamter aperture
	f4.g_convolvedflux_3_20_flux, f4.g_convolvedflux_3_20_fluxsigma,
	f4.r_convolvedflux_3_20_flux, f4.r_convolvedflux_3_20_fluxsigma,
	f4.i_convolvedflux_3_20_flux, f4.i_convolvedflux_3_20_fluxsigma,
	f4.z_convolvedflux_3_20_flux, f4.z_convolvedflux_3_20_fluxsigma,
	f4.y_convolvedflux_3_20_flux, f4.y_convolvedflux_3_20_fluxsigma,
	f4.g_convolvedflux_3_20_flag,
	f4.r_convolvedflux_3_20_flag,
	f4.i_convolvedflux_3_20_flag,
	f4.z_convolvedflux_3_20_flag,
	f4.y_convolvedflux_3_20_flag,
	-- PSF-corrected aperture photometry **before deblending**
	-- 2_15: 1.1 arcsec seeing; 1.5 arcsec diamter aperture
	f5.g_undeblended_convolvedflux_2_15_flux, f5.g_undeblended_convolvedflux_2_15_fluxsigma,
	f5.r_undeblended_convolvedflux_2_15_flux, f5.r_undeblended_convolvedflux_2_15_fluxsigma,
	f5.i_undeblended_convolvedflux_2_15_flux, f5.i_undeblended_convolvedflux_2_15_fluxsigma,
	f5.z_undeblended_convolvedflux_2_15_flux, f5.z_undeblended_convolvedflux_2_15_fluxsigma,
	f5.y_undeblended_convolvedflux_2_15_flux, f5.y_undeblended_convolvedflux_2_15_fluxsigma,
	f5.g_undeblended_convolvedflux_2_15_flag,
	f5.r_undeblended_convolvedflux_2_15_flag,
	f5.i_undeblended_convolvedflux_2_15_flag,
	f5.z_undeblended_convolvedflux_2_15_flag,
	f5.y_undeblended_convolvedflux_2_15_flag,
	-- 3_20: 1.3 arcsec seeing; 2.0 arcsec diamter aperture
	f5.g_undeblended_convolvedflux_3_20_flux, f5.g_undeblended_convolvedflux_3_20_fluxsigma,
	f5.r_undeblended_convolvedflux_3_20_flux, f5.r_undeblended_convolvedflux_3_20_fluxsigma,
	f5.i_undeblended_convolvedflux_3_20_flux, f5.i_undeblended_convolvedflux_3_20_fluxsigma,
	f5.z_undeblended_convolvedflux_3_20_flux, f5.z_undeblended_convolvedflux_3_20_fluxsigma,
	f5.y_undeblended_convolvedflux_3_20_flux, f5.y_undeblended_convolvedflux_3_20_fluxsigma,
	f5.g_undeblended_convolvedflux_3_20_flag,
	f5.r_undeblended_convolvedflux_3_20_flag,
	f5.i_undeblended_convolvedflux_3_20_flag,
	f5.z_undeblended_convolvedflux_3_20_flag,
	f5.y_undeblended_convolvedflux_3_20_flag,
	-- SDSS Shape without PSF correction (Using i-band; can use others too)
	f2.i_sdssshape_shape11 as i_sdss_shape_11,
	f2.i_sdssshape_shape12 as i_sdss_shape_12,
	f2.i_sdssshape_shape22 as i_sdss_shape_22,
	f2.i_sdssshape_shape11sigma as i_sdss_shape_11_err,
	f2.i_sdssshape_shape12sigma as i_sdss_shape_12_err,
	f2.i_sdssshape_shape22sigma as i_sdss_shape_22_err,
	-- Shape of the CModel model
	m.i_cmodel_exp_ellipse_11, m.i_cmodel_exp_ellipse_22, m.i_cmodel_exp_ellipse_12,
	m.i_cmodel_ellipse_11, m.i_cmodel_ellipse_22, m.i_cmodel_ellipse_12,
	m.r_cmodel_exp_ellipse_11, m.r_cmodel_exp_ellipse_22, m.r_cmodel_exp_ellipse_12,
	m.r_cmodel_ellipse_11, m.r_cmodel_ellipse_22, m.r_cmodel_ellipse_12,
	-- Flags for later selection
	-- 1. Source is outside usable exposure region
	f1.g_pixelflags_edge,
	f1.r_pixelflags_edge,
	f1.i_pixelflags_edge,
	f1.z_pixelflags_edge,
	f1.y_pixelflags_edge,
	-- 2. Saturated or interpolated pixels on the footprint (not center)
	f1.g_pixelflags_saturated,
	f1.r_pixelflags_saturated,
	f1.i_pixelflags_saturated,
	f1.z_pixelflags_saturated,
	f1.y_pixelflags_saturated,
	f1.g_pixelflags_interpolated,
	f1.r_pixelflags_interpolated,
	f1.i_pixelflags_interpolated,
	f1.z_pixelflags_interpolated,
	f1.y_pixelflags_interpolated,
	-- 3. Center is close to a clipped or suspicious pixel
	f1.g_pixelflags_suspectcenter,
	f1.r_pixelflags_suspectcenter,
	f1.i_pixelflags_suspectcenter,
	f1.z_pixelflags_suspectcenter,
	f1.y_pixelflags_suspectcenter,
	f1.g_pixelflags_clippedcenter,
	f1.r_pixelflags_clippedcenter,
	f1.i_pixelflags_clippedcenter,
	f1.z_pixelflags_clippedcenter,
	f1.y_pixelflags_clippedcenter,
	f1.g_pixelflags_bright_object,
	f1.r_pixelflags_bright_object,
	f1.i_pixelflags_bright_object,
	f1.z_pixelflags_bright_object,
	f1.y_pixelflags_bright_object,
	-- 4. Failed SDSS centroid algorithm
	-- It means the center is not well defined.
	f2.g_sdsscentroid_flag,
	f2.r_sdsscentroid_flag,
	f2.i_sdsscentroid_flag,
	f2.z_sdsscentroid_flag,
	f2.y_sdsscentroid_flag,
	-- Mizuki photo-z
	-- There are different estimators, can decide later
	pz.photoz_best, pz.photoz_median, pz.photoz_mean,
	-- Photo-z confidence
	pz.photoz_conf_best, pz.photoz_conf_median, pz.photoz_conf_mean,
	pz.photoz_err68_min, photoz_err68_max,
	pz.photoz_err95_min, photoz_err95_max,
	-- Photo-z risk
	pz.photoz_risk_best, pz.photoz_risk_median, pz.photoz_risk_mean,
	-- Mizuki also estimates a stellar mass from the template
	pz.stellar_mass, pz.stellar_mass_err68_min, pz.stellar_mass_err68_max

FROM
	pdr2_wide.forced as f1
	LEFT JOIN pdr2_wide.forced2 as f2 USING (object_id)
	LEFT JOIN pdr2_wide.forced4 as f4 USING (object_id)
	LEFT JOIN pdr2_wide.forced5 as f5 USING (object_id)
	LEFT JOIN pdr2_wide.meas as m USING (object_id)
	LEFT JOIN pdr2_wide.photoz_mizuki as pz USING (object_id)

WHERE
    f1.ra between 149.38 and 150.8 
	AND f1.dec between 1.6 and 2.83
	AND f1.isprimary
	AND f1.r_cmodel_mag <= 23.6
	AND f1.g_cmodel_mag - f1.i_cmodel_mag BETWEEN 0 AND 1.6
	AND f1.r_cmodel_mag - f1.y_cmodel_mag BETWEEN -0.5 AND 1.2
	AND f1.g_inputcount_value >= 2
	AND f1.r_inputcount_value >= 2
	AND f1.i_inputcount_value >= 3
	AND f1.z_inputcount_value >= 3
	AND f1.y_inputcount_value >= 3
	AND f1.r_extendedness_value > 0
	AND f1.i_extendedness_value > 0
	AND NOT f1.g_pixelflags_saturatedcenter
	AND NOT f1.r_pixelflags_saturatedcenter
	AND NOT f1.i_pixelflags_saturatedcenter
	AND NOT f1.z_pixelflags_saturatedcenter
	AND NOT f1.y_pixelflags_saturatedcenter
	AND NOT f1.g_pixelflags_interpolatedcenter
	AND NOT f1.r_pixelflags_interpolatedcenter
	AND NOT f1.i_pixelflags_interpolatedcenter
	AND NOT f1.z_pixelflags_interpolatedcenter
	AND NOT f1.y_pixelflags_interpolatedcenter
	AND NOT f1.g_pixelflags_bad
	AND NOT f1.r_pixelflags_bad
	AND NOT f1.i_pixelflags_bad
	AND NOT f1.z_pixelflags_bad
	AND NOT f1.y_pixelflags_bad