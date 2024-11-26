
# -*- coding: utf-8 -*-
"""
Author:  Magnier Morgane (morgane.magnier@vattenfall.com) 

Copyright © 2024 Magnier Morgane
This work is part of a thesis project. The copyright of the thesis itself belongs to the student Magnier Morgane. 
Vattenfall has the right to use the thesis work, including its findings, methods, and conclusions, in its operations. Any material generated within the framework of this thesis that is subject to intellectual property protection (e.g., source code, computer program, design, or invention) is the property of Vattenfall, unless otherwise agreed in writing. Permission is hereby granted to view, copy, and share this work for educational or personal purposes only, provided that this copyright notice and this permission notice appear in all copies.
"""

import geemap
import ee

def diff_image(roi, band, image_1, image_2):

    band_img_1 = image_1.select(band)
    band_img_2 = image_2.select(band)

    diff = band_img_2.subtract(band_img_1).rename('diff')

    max = ee.Number(diff.reduceRegion(reducer = ee.Reducer.stdDev(), geometry = roi, scale = 10).get('diff'))
    min = max.multiply(-1)

    return diff

def detection_changes(roi, diff,thresh): 

    #thresholdGain = ee.Number(diff.reduceRegion(reducer = ee.Reducer.stdDev(), geometry = roi, scale = 10).get('diff'))
    thresholdGain = thresh
    thresholdLoss = -thresh

    diffClassified = ee.Image(0)

    diffClassified = diffClassified.where(diff.lte(thresholdLoss), 2)
    diffClassified = diffClassified.where(diff.gte(thresholdGain), 1)

    return diffClassified

def getChanges(roi, thresh, band, image_1, image_2): 

    diff = diff_image(roi, band, image_1, image_2)

    diffClassified = detection_changes(roi, diff,thresh)

    return diff, diffClassified
