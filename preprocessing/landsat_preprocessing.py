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

def mask_landsat_clouds(image):
  """Masks clouds in a andsat image using the QA band.

  Args:
      image (ee.Image): A Landsat image.

  Returns:
      ee.Image: A cloud-masked Landsat image.
  """
  qa = image.select('QA_PIXEL')

  # Bits 10 and 11 are clouds and cirrus, respectively.
  cloud_bit_mask = 1 << 2
  cirrus_bit_mask = 1 << 3
  cloud_shadow_bit_mask = 1 << 4

  # Both flags should be set to zero, indicating clear conditions.
  mask = (
      qa.bitwiseAnd(cloud_bit_mask)
      .eq(0)
      .And(qa.bitwiseAnd(cirrus_bit_mask).eq(0).And(qa.bitwiseAnd(cloud_shadow_bit_mask).eq(0)))
  )

  return image.updateMask(mask)

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

def get_l9_cloud_free_col(roi, indices, thresh):

    l9 = ee.ImageCollection("LANDSAT/LC09/C02/T1_L2").filterBounds(roi).filter(ee.Filter.lessThan('CLOUD_COVER', 80))
    l9 = l9.spectralIndices(indices)
    l9_masked = l9.map(mask_landsat_clouds)
    l9 = filter_col(l9_masked, roi, 'SR_B1',thresh)

    return l9


def get_l9_cloud_free_col_dates(roi, indices, thresh, start_date, end_date):

    l9 = ee.ImageCollection("LANDSAT/LC09/C02/T1_L2").filterBounds(roi).filter(ee.Filter.lessThan('CLOUD_COVER', 80)).filterDate(start_date, end_date)
    l9 = l9.spectralIndices(indices)
    l9_masked = l9.map(mask_landsat_clouds)
    l9 = filter_col(l9_masked, roi, 'SR_B1', thresh)

    return l9


def get_l8_cloud_free_col(roi, indices, thresh):

    l8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2").filterBounds(roi).filter(ee.Filter.lessThan('CLOUD_COVER', 80))
    l8 = l8.spectralIndices(indices)
    l8_masked = l8.map(mask_landsat_clouds)
    l8 = filter_col(l8_masked, roi, 'SR_B1',thresh)

    return l8


def get_l8_cloud_free_col_dates(roi, indices, thresh, start_date, end_date):

    l8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2").filterBounds(roi).filter(ee.Filter.lessThan('CLOUD_COVER', 80)).filterDate(start_date, end_date)
    l8 = l8.spectralIndices(indices)
    l8_masked = l8.map(mask_landsat_clouds)
    l8 = filter_col(l8_masked, roi, 'SR_B1',thresh)

    return l8


def get_l7_cloud_free_col(roi, indices, thresh):

    l7 = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2").filterBounds(roi).filter(ee.Filter.lessThan('CLOUD_COVER', 80))
    l7 = l7.spectralIndices(indices)
    l7_masked = l7.map(mask_landsat_clouds)
    l7 = filter_col(l7_masked, roi, 'SR_B1',thresh)

    return l7


def get_l7_cloud_free_col_dates(roi, indices, thresh, start_date, end_date):

    l7 = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2").filterBounds(roi).filter(ee.Filter.lessThan('CLOUD_COVER', 80)).filterDate(start_date, end_date)
    l7 = l7.spectralIndices(indices)
    l7_masked = l7.map(mask_landsat_clouds)
    l7 = filter_col(l7_masked, roi, 'SR_B1', thresh)

    return l7


def get_l5_cloud_free_col(roi, indices, thresh):

    l5 = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2").filterBounds(roi).filter(ee.Filter.lessThan('CLOUD_COVER', 80))
    l5 = l5.spectralIndices(indices)
    l5_masked = l5.map(mask_landsat_clouds)
    l5 = filter_col(l5_masked, roi, 'SR_B1',thresh)

    return l5


def get_l5_cloud_free_col_dates(roi, indices, thresh, start_date, end_date):

    l5 = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2").filterBounds(roi).filter(ee.Filter.lessThan('CLOUD_COVER', 80)).filterDate(start_date, end_date)
    l5 = l5.spectralIndices(indices)
    l5_masked = l5.map(mask_landsat_clouds)
    l5 = filter_col(l5_masked, roi, 'SR_B1', thresh)

    return l5



'''
 def get_l4_cloud_free_col(roi, indices, thresh):

   l4 = ee.ImageCollection("LANDSAT/LT04/C02/T1_L2").filterBounds(roi).filter(ee.Filter.lessThan('CLOUD_COVER', 80))
   l4 = l.spectralIndices(indices)
   l4_masked = l4.map(mask_landsat_clouds)
   l4 = filter_col(l4_masked, roi, 'SR_B1',thresh)

   return l4
'''

 
