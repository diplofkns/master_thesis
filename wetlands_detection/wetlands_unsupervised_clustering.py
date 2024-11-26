# -*- coding: utf-8 -*- 

""" 
Author:  Magnier Morgane (morgane.magnier@vattenfall.com) 

Copyright © 2024 Magnier Morgane
This work is part of a thesis project. The copyright of the thesis itself belongs to the student Magnier Morgane. 
Vattenfall has the right to use the thesis work, including its findings, methods, and conclusions, in its operations. Any material generated within the framework of this thesis that is subject to intellectual property protection (e.g., source code, computer program, design, or invention) is the property of Vattenfall, unless otherwise agreed in writing. Permission is hereby granted to view, copy, and share this work for educational or personal purposes only, provided that this copyright notice and this permission notice appear in all copies.
"""


import ee
import geemap
import eemont

ee.Authenticate()
ee.Initialize()


def clustering(image,roi): 

    bands = bands = [
    "ANDWI",
    "AWEInsh",
    "AWEIsh",
    "LSWI",
    "MBWI",
    "MLSWI26",
    "MLSWI27",
    "MNDWI",
    "MuWIR",
    "NDCI",
    "NDPonI",
    "NDTI",
    "NDVIMNDWI",
    "NDWI",
    "NDWIns",
    "NWI",
    "S2WI",
    "SWM",
    "WI1",
    "WI2",
    "WI2015",
    "WRI"
]
    
    training_samples = image.select(bands).sample(region=roi, scale=5, numPixels=15000, geometries=True)
    
    clusterer = ee.Clusterer.wekaKMeans(nClusters = 2, distanceFunction="Euclidean", maxIterations=75).train(training_samples)

    image_ndwi_clusters = image.addBands(image.cluster(clusterer)).select(['NDWI', 'cluster'])

    cluster = image.addBands(image.cluster(clusterer).rename('cluster'))

    return cluster

def water_classification(image, roi):
    
    image = clustering(image, roi)

    cluster_1 = image.updateMask(image.select('cluster').eq(1))
    cluster_2 = image.updateMask(image.select('cluster').eq(0))
    
    mean_1 = cluster_1.select('B1').reduceRegion(ee.Reducer.mean(),geometry = image.geometry(), scale = 10)
    mean_2 = cluster_2.select('B1').reduceRegion(ee.Reducer.mean(),geometry = image.geometry(), scale = 10)
    
    mean_1_value = mean_1.getNumber('B1')
    mean_2_value = mean_2.getNumber('B1')
    
    #display(mean_1_value)
    #display(mean_2_value)
    
    water = ee.Image(ee.Algorithms.If(mean_1_value.gt(mean_2_value), image.select('cluster').eq(0), image.select('cluster').eq(1)))

    return water

def apply_mask(image,wetlands_mask):
    image = image.updateMask(wetlands_mask)
    return image

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

def getWetlandsS2(roi, min_water_date, max_water_date):
    
    s2_col = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED').filterBounds(roi).spectralIndices('water')

    flooded_img = filter_col(s2_col.filterDate(max_water_date,max_water_date.advance(3, 'day')),roi, 'B1', 0.95).first().resample('bicubic').clip(roi)
    dry_img = filter_col(s2_col.filterDate(min_water_date, min_water_date.advance(3, 'day')),roi, 'B1', 0.95).first().resample('bicubic').clip(roi)
    display(flooded_img)
    display(dry_img)
    #flooded_img = get_s2_clouds_free_col(roi, max_water_date, max_water_date.advance(1,'day')).spectralIndices('water').first().resample('bicubic').clip(roi)
    #dry_img = get_s2_clouds_free_col(roi, min_water_date, min_water_date.advance(1,'day')).spectralIndices('water').first().resample('bicubic').clip(roi)

    flooded_mask = water_classification(flooded_img, roi)
    dry_mask = water_classification(dry_img, roi)

    wetlands_mask = flooded_mask.subtract(dry_mask)
    
    return wetlands_mask
