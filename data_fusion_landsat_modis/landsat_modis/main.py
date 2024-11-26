# -*- coding: utf-8 -*-
"""
Author: Morgane Magnier (morgane.magnier@vattenfall.com)

Copyright © [2024] [Magnier Morgane]  

This work is part of a thesis project. The copyright of the thesis itself belongs to the student Morgane Magnier. 

Vattenfall has the right to use the thesis work, including its findings, methods, and conclusions, in its operations. 
Any material generated within the framework of this thesis that is subject to intellectual property protection 
(e.g., source code, computer program, design, or invention) is the property of Vattenfall, unless otherwise agreed in writing.

Permission is hereby granted to view, copy, and share this work for educational or personal purposes only, provided that this
copyright notice and this permission notice appear in all copies.
"""


### GLOBALS ##########################################################################################################

import ee, geemap, eemont
from GEE_ImageFusion import *

ee.Initialize()
ee.Authenticate()


def data_fusion(startDate, endDate): 
    # region of interest
    # location for prediction (outside scene overlap areas)
    region = ee.Geometry.Polygon([[[17.204933,60.402663],[17.204933,60.455525],[17.2645,60.455525],[17.2645,60.402663],[17.204933,60.402663]]])

    # common bands between sensors that will be used for fusion
    #   would need to add functions for other indices (evi etc.).
    #   ndvi is exported as 16 bit int so would have to make sure not to rescale
    #   reflectance bands if that is what you are planning to predict (line 206)
    #   some resturcturing of this code would be necessary if the goal was to
    #   predict a multiband image but the core functions should work for any number
    #   of bands
    commonBandNames = ee.List(['ndvi'])

    # image collections to use in fusion
    # NOTE: if using older Landsat and not using NDVI one would have to modify the
    # get_paired_collections script because this script harmonizes NDVI from
    # L5 & L7 to L8 based on Roy et al. 2016 (see etmToOli and getPaired functions)
    landsatCollection = 'LANDSAT/LC08/C02/T1_L2'
    modisCollection = 'MODIS/006/MCD43A4'

    # landsat band names including qc band for masking
    bandNamesLandsat = ee.List(['blue', 'green', 'red',
                                'nir', 'swir1', 'swir2', 'pixel_qa'])
    landsatBands = ee.List([1, 2, 3, 4, 5, 6, 10])

    # modis band names
    bandNamesModis = ee.List(['blue', 'green', 'red', 'nir', 'swir1', 'swir2'])
    modisBands = ee.List([2, 3, 0, 1, 5, 6])

    # radius of moving window
    # Note: Generally, larger windows are better but as the window size increases,
    # so does the memory requirement and we quickly will surpass the memory
    # capacity of a single node (in testing 13 was max size for single band, and
    # 10 was max size for up to 6 bands)
    kernelRadius = ee.Number(10)
    kernel = ee.Kernel.square(kernelRadius)
    numPixels = kernelRadius.add(kernelRadius.add(1)).pow(2)

    # number of land cover classes in scene
    coverClasses = 7
    # to export the images to an asset we need the path to the assets folder
    path = 'projects/ee-magniermorgane/assets/'
    scene_name = 'NDVI_'

    ### GET FILTERED COL ##########################################################################################################

    # sorted, filtered, paired image retrieval
    paired = getPaired(startDate, endDate,
                    landsatCollection, landsatBands, bandNamesLandsat,
                    modisCollection, modisBands, bandNamesModis,
                    commonBandNames, region)

    subs = makeSubcollections(paired)
    # subs_meta = subs.getInfo()


    ### PREDICT ##########################################################################################################

    # loop through each list of paired images
    num_lists = subs.length().getInfo()
    for i in range(0, num_lists):
        # determine the number of modis images between pairs
        num_imgs = ee.List(ee.List(subs.get(i)).get(2)).length()

        # determine the remainder of images, if not groups of 10
        remaining = num_imgs.mod(10)

        # create sequence of starting indices for modis images
        index_seq = ee.List.sequence(0, num_imgs.subtract(remaining), 10)

        # images to be grouped and predicted
        subList = ee.List(ee.List(subs.get(i)).get(2))

        # loop through indices predicting in batches of 10
        for x in range(0, index_seq.length().getInfo()):
            # starting index
            start = ee.Number(index_seq.get(x))

            # ending index
            end = ee.Algorithms.If(start.add(10).gt(num_imgs),
                                num_imgs,
                                start.add(10))

            # group of images to predict
            pred_group = subList.slice(start, end)
            landsat_t01 = ee.List(ee.List(subs.get(i)).get(0))
            modis_t01 = ee.List(ee.List(subs.get(i)).get(1))
            modis_tp = pred_group

            # get the start and end day values and year for the group to use
            # to label the file when exported to asset
            startDay = ee.Number.parse(ee.ImageCollection(pred_group)
                                    .first()
                                    .get('DOY'))
            endDay = ee.Number.parse(ee.ImageCollection(pred_group)
                                    .sort('system:time_start', False)
                                    .first()
                                    .get('DOY'))
            year = ee.Date(ee.ImageCollection(pred_group)
                        .sort('system:time_start', False)
                        .first()
                        .get('system:time_start')).format('Y')

            # start and end day of year
            doys = landsat_t01 \
                .map(lambda img: ee.String(ee.Image(img).get('DOY')).cat('_'))

            # register images
            landsat_t01, modis_t01, modis_tp = registerImages(landsat_t01,
                                                            modis_t01,
                                                            modis_tp)

            # prep landsat imagery (mask and format)
            maskedLandsat, pixPositions, pixBN = prepLandsat(landsat_t01,
                                                            kernel,
                                                            numPixels,
                                                            commonBandNames,
                                                            doys,
                                                            coverClasses)

            # prep modis imagery (mask and format)
            modSorted_t01, modSorted_tp = prepMODIS(modis_t01, modis_tp, kernel,
                                                    numPixels, commonBandNames,
                                                    pixBN)

            # calculate spectral distance
            specDist = calcSpecDist(maskedLandsat, modSorted_t01,
                                    numPixels, pixPositions)

            # calculate spatial distance
            spatDist = calcSpatDist(pixPositions)

            # calculate weights from the spatial and spectral distances
            weights = calcWeight(spatDist, specDist)

            # calculate the conversion coefficients
            coeffs = calcConversionCoeff(maskedLandsat, modSorted_t01,
                                        doys, numPixels, commonBandNames)

            # predict all modis images in modis tp collection
            prediction = modSorted_tp \
                .map(lambda image:
                    predictLandsat(landsat_t01, modSorted_t01,
                                    doys, ee.List(image),
                                    weights, coeffs,
                                    commonBandNames, numPixels))

            # create a list of new band names to apply to the multiband ndvi image
            # NOTE: cant export with names starting with 0
            preds = ee.ImageCollection(prediction).toBands()
            dates = modis_tp.map(lambda img:
                                ee.Image(img).get('system:time_start'))
            predNames = ee.List.sequence(0, prediction.length().subtract(1)) \
                .map(lambda i:
                    commonBandNames\
                        .map(lambda name:
                            ee.String(name)
                            .cat(ee.String(ee.Number(dates.get(i)).format()))))\
                .flatten()

            # export all predictions as a single multiband image
            # each band name corresponds to the timestamp for the image
            task = ee.batch.Export.image.toAsset(
                image=preds.rename(predNames).multiply(10000).toInt16(),
                description=ee.String(scene_name)
                            .cat(year)
                            .cat('_')
                            .cat(startDay.format())
                            .cat('_').cat(endDay.format()).getInfo(),
                assetId=ee.String(path)
                        .cat(ee.String(scene_name))
                        .cat(year)
                        .cat('_')
                        .cat(startDay.format())
                        .cat('_')
                        .cat(endDay.format()).getInfo(),
                region=ee.Image(prediction.get(0)).geometry(),
                scale=30)
            
            return task 
            
