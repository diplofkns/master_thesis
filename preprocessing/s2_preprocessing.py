# -*- coding: utf-8 -*-

""" 
Author:  Magnier Morgane (morgane.magnier@vattenfall.com) 

Copyright © 2024 Magnier Morgane
This work is part of a thesis project. The copyright of the thesis itself belongs to the student Magnier Morgane. 
Vattenfall has the right to use the thesis work, including its findings, methods, and conclusions, in its operations.
Any material generated within the framework of this thesis that is subject to intellectual property protection (e.g., source code, computer program, design, or invention) is the property of Vattenfall, unless otherwise agreed in writing. Permission is hereby granted to view, copy, and share this work for educational or personal purposes only, provided that this copyright notice and this permission notice appear in all copies.
"""

import ee,eemont

ee.Authenticate()
ee.Initialize()

def get_s2_sr_cld_col(roi):
    # Import and filter S2 SR.
    CLOUD_FILTER = 80 
    s2_sr_col = (ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
        .filterBounds(roi)
        .filter(ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE', CLOUD_FILTER)))

    # Import and filter s2cloudless.
    s2_cloudless_col = ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY').filterBounds(roi)

    # Join the filtered s2cloudless collection to the SR collection by the 'system:index' property.
    return ee.ImageCollection(ee.Join.saveFirst('s2cloudless').apply(**{
        'primary': s2_sr_col,
        'secondary': s2_cloudless_col,
        'condition': ee.Filter.equals(**{
            'leftField': 'system:index',
            'rightField': 'system:index'
        })
    }))

def add_cloud_bands(img):
    CLD_PRB_THRESH = 40

    # Get s2cloudless image, subset the probability band.
    cld_prb = ee.Image(img.get('s2cloudless')).select('probability')

    # Condition s2cloudless by the probability threshold value.
    is_cloud = cld_prb.gt(CLD_PRB_THRESH).rename('clouds')

    # Add the cloud probability layer and cloud mask as image bands.
    return img.addBands(ee.Image([cld_prb, is_cloud]))

def add_shadow_bands(img):

    NIR_DRK_THRESH = 0.15
    CLD_PRJ_DIST = 2

    # Identify water pixels from the SCL band.
    not_water = img.select('SCL').neq(6)

    # Identify dark NIR pixels that are not water (potential cloud shadow pixels).
    SR_BAND_SCALE = 1e4
    dark_pixels = img.select('B8').lt(NIR_DRK_THRESH*SR_BAND_SCALE).multiply(not_water).rename('dark_pixels')

    # Determine the direction to project cloud shadow from clouds (assumes UTM projection).
    shadow_azimuth = ee.Number(90).subtract(ee.Number(img.get('MEAN_SOLAR_AZIMUTH_ANGLE')))

    # Project shadows from clouds for the distance specified by the CLD_PRJ_DIST input.
    cld_proj = (img.select('clouds').directionalDistanceTransform(shadow_azimuth, CLD_PRJ_DIST*10)
        .reproject(**{'crs': img.select(0).projection(), 'scale': 100})
        .select('distance')
        .mask()
        .rename('cloud_transform'))

    # Identify the intersection of dark pixels with cloud shadow projection.
    shadows = cld_proj.multiply(dark_pixels).rename('shadows')

    # Add dark pixels, cloud projection, and identified shadows as image bands.
    return img.addBands(ee.Image([dark_pixels, cld_proj, shadows]))

def add_cld_shdw_mask(img):
    BUFFER = 100
    # Add cloud component bands.
    img_cloud = add_cloud_bands(img)

    # Add cloud shadow component bands.
    img_cloud_shadow = add_shadow_bands(img_cloud)

    # Combine cloud and shadow mask, set cloud and shadow as value 1, else 0.
    is_cld_shdw = img_cloud_shadow.select('clouds').add(img_cloud_shadow.select('shadows')).gt(0)

    # Remove small cloud-shadow patches and dilate remaining pixels by BUFFER input.
    # 20 m scale is for speed, and assumes clouds don't require 10 m precision.
    is_cld_shdw = (is_cld_shdw.focalMin(2).focalMax(BUFFER*2/20)
        .reproject(**{'crs': img.select([0]).projection(), 'scale': 20})
        .rename('cloudmask'))

    # Add the final cloud-shadow mask to the image.
    return img_cloud_shadow.addBands(is_cld_shdw)

def apply_cld_shdw_mask(img):
    # Subset the cloudmask band and invert it so clouds/shadow are 0, else 1.
    not_cld_shdw = img.select('cloudmask').Not()

    # Subset reflectance bands and update their masks, return the result.
    return img.select('B.*').updateMask(not_cld_shdw)

def filter_col(col, roi, band, thresh):
    
    col = col.map(lambda image : image.clip(roi))

    def count_pixels(image,roi): 
        pixel_count = image.select(band).reduceRegion(
            reducer= ee.Reducer.count(),
            geometry=roi,
            scale=10,
            maxPixels=1e9
        ).get(band)
        return image.set('pixel_count', pixel_count)

    nb_pixels_ts = col.map(lambda image: count_pixels(image, roi))

    # Get the image with the maximum pixel count
    max_pixel_count_image = nb_pixels_ts.sort('pixel_count', False).first()
    ref_img_pixel_count = max_pixel_count_image.get('pixel_count').getInfo()
    pixel_count_threshold = ref_img_pixel_count * thresh

    # Filter the collection based on the pixel count threshold
    filtered_col = nb_pixels_ts.filter(ee.Filter.gte('pixel_count', pixel_count_threshold))

    return filtered_col

def get_s2_cloud_free_col(roi, thresh, indices = None):

    s2_sr = get_s2_sr_cld_col(roi)
    s2_sr_masked = (s2_sr.map(add_cld_shdw_mask)
                                .map(apply_cld_shdw_mask))

    s2 = filter_col(s2_sr_masked, roi, 'B1', thresh)

    if indices : 
        s2 = s2.spectralIndices(indices)
        
    return s2

def get_s2_cloud_free_col_dates(roi, thresh, start_date,end_date,indices = None):

    s2_sr = get_s2_sr_cld_col(roi).filterDate(start_date,end_date)
    s2_sr_masked = (s2_sr.map(add_cld_shdw_mask)
                                .map(apply_cld_shdw_mask))

    s2 = filter_col(s2_sr_masked, roi, 'B1', thresh)

    if indices : 
        s2 = s2.spectralIndices(indices)

    return s2
