
# -*- coding: utf-8 -*-

""" 
Author:  Magnier Morgane (morgane.magnier@vattenfall.com) 

Copyright © 2024 Magnier Morgane
This work is part of a thesis project. The copyright of the thesis itself belongs to the student Magnier Morgane. 
Vattenfall has the right to use the thesis work, including its findings, methods, and conclusions, in its operations.
Any material generated within the framework of this thesis that is subject to intellectual property protection (e.g., source code, computer program, design, or invention) is the property of Vattenfall, unless otherwise agreed in writing. Permission is hereby granted to view, copy, and share this work for educational or personal purposes only, provided that this copyright notice and this permission notice appear in all copies.
"""

import ee, eemont

ee.Authenticate()
ee.Initialize()

def maskMODIS(image):
    """
    Mask snow covered and extremely high albedo areas from the modis images.

    Parameters
    ----------
    image : image.Image
        MODIS image.

    Returns
    -------
    image.image
        Masked MODIS image.

    """
    # calculate snow water index for the image
    swi = image.expression(
        '(green * (nir - swir1)) / ((green + nir) * (nir + swir1))',
        {'green': image.select(['green']),
         'nir': image.select(['nir']),
         'swir1': image.select(['swir1'])
         }).rename('swi')

    # mask out values of swi above 0.1
    mask = swi.lt(0.1)

    return image \
        .updateMask(mask) \
        .copyProperties(image, ['system:time_start', 'system:id'])

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


def get_modis_cloud_free_col(roi, thresh, indices = None):

    modis_raw = ee.ImageCollection('MODIS/006/MCD43A4').filterBounds(roi)
    modis_masked = modis_raw.map(maskMODIS)

    modis = filter_col(modis_masked, roi, 'Nadir_Reflectance_Band1', thresh)

    if indices : 
        modis = modis.spectralIndices(indices)

    return modis

def get_modis_cloud_free_col_dates(roi, thresh, start_date,end_date, indices = None):

    bandNamesModis = ee.List(['blue', 'green', 'red', 'nir', 'swir1', 'swir2'])
    modisBands = ee.List([2, 3, 0, 1, 5, 6])
    modis_raw = ee.ImageCollection('MODIS/006/MCD43A4') \
            .filterBounds(roi)\
            .select(modisBands, bandNamesModis) \
            .map(maskMODIS) \
              
    modis_masked = modis_raw.map(maskMODIS)

    modis = filter_col(modis_masked, roi, 'Nadir_Reflectance_Band1', thresh)

    if indices : 
        modis = modis.spectralIndices(indices)

    return modis
