# -*- coding: utf-8 -*- 

""" 
Author:  Magnier Morgane (morgane.magnier@vattenfall.com) 

Copyright © 2024 Magnier Morgane
This work is part of a thesis project. The copyright of the thesis itself belongs to the student Magnier Morgane. 
Vattenfall has the right to use the thesis work, including its findings, methods, and conclusions, in its operations. Any material generated within the framework of this thesis that is subject to intellectual property protection (e.g., source code, computer program, design, or invention) is the property of Vattenfall, unless otherwise agreed in writing. Permission is hereby granted to view, copy, and share this work for educational or personal purposes only, provided that this copyright notice and this permission notice appear in all copies.
"""

import ee, eemont, geemap
import sys 
sys.path.append('../utils')
import l8_selecting_clouds_free

def add_elev_slope(image): 

    jax_dsm = ee.ImageCollection('JAXA/ALOS/AW3D30/V3_2')
    jax_elevation = jax_dsm.select('DSM')

    proj = jax_elevation.first().select(0).projection()
    slopeReprojected =  (jax_elevation.mosaic() \
                                .setDefaultProjection(proj)).resample('bicubic')

    # Reduce the collection with a median reducer.
    elevation = slopeReprojected.reduce(ee.Reducer.mean()).rename('elev')

    slope = ee.Terrain.slope(elevation).rename('slope')

    image = image.addBands(elevation).addBands(slope)

    return image

def normalize(image, roi):
    bandNames = image.bandNames()
  # Compute min and max of the image
    minDict = image.reduceRegion(
        reducer = ee.Reducer.min(),
        geometry = roi,
        scale = 10,
        maxPixels=1e9,
        bestEffort=True,
        tileScale= 16)
    
    maxDict = image.reduceRegion(
        reducer= ee.Reducer.max(),
        geometry= roi,
        scale= 10,
        maxPixels= 1e9,
        bestEffort= True,
        tileScale= 16)
    
    mins = ee.Image.constant(minDict.values(bandNames))
    maxs = ee.Image.constant(maxDict.values(bandNames))

    normalized = image.subtract(mins).divide(maxs.subtract(mins))
    return normalized

def create_composite_per_month(roi, indices, year, month, resample=True):

    start_date = ee.Date.fromYMD(year,month,1)
    end_date = start_date.advance(1, 'month')

    l8_col = l8_selecting_clouds_free.get_clouds_free_collection_in_roi(roi,start_date, end_date)

    l8_col_index = l8_col.spectralIndices(indices)

    if resample:
        l8_rs = l8_col_index.map(lambda i: i.resample('bicubic'))

    img = l8_rs.median()

    composite = add_elev_slope(img)

    composite = normalize(composite, roi)

    return composite.set({
        'year': year,
        'month': month,
        'date': ee.Date.fromYMD(year, month, 1)
    })

def create_composite_per_year(roi, indices, year, resample=True): 

    start_date = ee.Date.fromYMD(year, 5 ,1)
    end_date = start_date.advance(4, 'month')

    l8_col = l8_selecting_clouds_free.get_clouds_free_collection_in_roi(roi,start_date, end_date)

    l8_col_index = l8_col.spectralIndices(indices)

    if resample:
        l8_rs = l8_col_index.map(lambda i: i.resample('bicubic'))

    img = l8_rs.median()

    composite = add_elev_slope(img)

    composite = normalize(composite, roi)

    return composite.set({
        'year': year,
    })


def create_composite_per_date(roi, indices, start_date, end_date, resample=True): 

    l8_col = l8_selecting_clouds_free.get_clouds_free_collection_in_roi(roi,start_date, end_date)

    l8_col_index = l8_col.spectralIndices(indices)

    if resample:
        l8_rs = l8_col_index.map(lambda i: i.resample('bicubic'))

    img = l8_rs.median()

    composite = add_elev_slope(img)

    composite = normalize(composite, roi)

    return composite.set({
        'year': year,
        'month': month,
        'date': ee.Date.fromYMD(year, month, 1)
    })


def create_composite_single_image(image,roi):

    img = image.resample('bicubic')

    composite = add_elev_slope(img)

    composite = normalize(composite, roi)

    return composite
