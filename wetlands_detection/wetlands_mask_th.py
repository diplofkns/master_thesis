# -*- coding: utf-8 -*- 

""" 
Author:  Magnier Morgane (morgane.magnier@vattenfall.com) 

Copyright © 2024 Magnier Morgane
This work is part of a thesis project. The copyright of the thesis itself belongs to the student Magnier Morgane. 
Vattenfall has the right to use the thesis work, including its findings, methods, and conclusions, in its operations. Any material generated within the framework of this thesis that is subject to intellectual property protection (e.g., source code, computer program, design, or invention) is the property of Vattenfall, unless otherwise agreed in writing. Permission is hereby granted to view, copy, and share this work for educational or personal purposes only, provided that this copyright notice and this permission notice appear in all copies.
"""

import ee

def apply_ndwi(img):
    exp = '(G-NIR)/(G+NIR)'
    ndwi = (ee.Image().expression(expression=exp,opt_map={'NIR': img.select('B5'), 'G': img.select('B3')},)).rename('NDWI')
    return img.addBands(ndwi)  

def water_classification_th(image,th):
    water = image.select('NDWI').gt(th).rename('water_mask')
    return water

def apply_mask(image, col, min_water_date, max_water_date, th = -0.15):
    wetlands_mask = getWetlands(col, min_water_date, max_water_date, th)
    image = image.updateMask(wetlands_mask)
    return image

def getWetlands(col, min_water_date, max_water_date, th): 

    flooded_img = col.filterDate(max_water_date, max_water_date.advance(1, 'day')).first()
    dry_img = col.filterDate(min_water_date, min_water_date.advance(1, 'day')).first()

    flooded_img = apply_ndwi(flooded_img)
    dry_img = apply_ndwi(dry_img)

    flooded_mask = water_classification_th(flooded_img,th)
    dry_mask = water_classification_th(dry_img,th)

    wetlands_mask = flooded_mask.subtract(dry_mask)

    return wetlands_mask
