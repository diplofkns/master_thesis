
# -*- coding: utf-8 -*- 

""" 
Author:  Magnier Morgane (morgane.magnier@vattenfall.com) 

Copyright © 2024 Magnier Morgane
This work is part of a thesis project. The copyright of the thesis itself belongs to the student Magnier Morgane. 
Vattenfall has the right to use the thesis work, including its findings, methods, and conclusions, in its operations. Any material generated within the framework of this thesis that is subject to intellectual property protection (e.g., source code, computer program, design, or invention) is the property of Vattenfall, unless otherwise agreed in writing. Permission is hereby granted to view, copy, and share this work for educational or personal purposes only, provided that this copyright notice and this permission notice appear in all copies.
"""

import ee, eemont
import geemap
import math
import pandas as pd

import matplotlib.pyplot as plt
ee.Authenticate()
ee.Initialize()

def water_cluster(s1, start_date, end_date, wind_filter):

    ##### GENERAL - TO BE CHANGE IF NEEDED #################################################################################
    roi =  ee.Geometry.Polygon([[[17.204933,60.402663],[17.204933,60.455525],[17.2645,60.455525],[17.2645,60.402663],[17.204933,60.402663]]])

    # Filter the Sentinel-1 collection by metadata properties.
    vv_vh_iw = (
        s1
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
        .filter(ee.Filter.eq('instrumentMode', 'IW'))
    )

    # Separate ascending and descending orbit images into distinct collections.
    vv_vh_iw_asc = vv_vh_iw.filter(
        ee.Filter.eq('orbitProperties_pass', 'ASCENDING')
    )
    vv_vh_iw_desc = vv_vh_iw.filter(
        ee.Filter.eq('orbitProperties_pass', 'DESCENDING')
    )

    ##### WIND FILTER #################################################################################

    if wind_filter == True : 

        def create_wind_band(image_collection):
            # Get date property from the image collection, format it
            def format_date(image):
                return image.set('date', image.date().format('Y-MM-dd'))

            dates = image_collection.map(format_date)

            # Get a list of the dates.
            dates_list = dates.aggregate_array('date')

            # Function to get wind data for each date
            def get_wind_data(date):
                wx = ee.ImageCollection('NOAA/CFSV2/FOR6H').filterDate(ee.Date(date), ee.Date(date).advance(1, 'day'))
                v_wind = wx.select(['v-component_of_wind_height_above_ground']).max()
                u_wind = wx.select(['u-component_of_wind_height_above_ground']).max()
                
                wind_speed = v_wind.pow(2).add(u_wind.pow(2)).sqrt().multiply(3.6)  # Convert from m/s to km/h
                return wind_speed.rename('windy').set('date', date)
            
            # Map the wind data function over the dates list
            wind_speed_images = ee.ImageCollection(dates_list.map(get_wind_data))

            # Define an inner join.
            inner_join = ee.Join.inner()

            # Specify an equals filter for image dates
            filter_date_eq = ee.Filter.equals(leftField='date', rightField='date')

            # Apply the join to combine wind data with the image collection
            inner_joined = inner_join.apply(dates, wind_speed_images, filter_date_eq)

            # Function to concatenate the joined images
            def concatenate_images(joined):
                return ee.Image.cat(joined.get('primary'), joined.get('secondary'))

            # Map the concatenate function over the joined collection
            joined_collection = ee.ImageCollection(inner_joined.map(concatenate_images))

            return joined_collection

        def filter_images_by_wind_speed(image_collection, roi , wind_speed_threshold=12.0, min_pixels=1000):
            # Create the wind band and join it to the image collection
            image_collection_with_wind = create_wind_band(image_collection)

            # Apply the wind mask to the image collection
            def wind_mask(image):
                # Get the wind speed from the image
                wind_speed_image = image.select('windy')
                
                # Create a mask: 1 where wind speed < threshold, 0 otherwise
                mask = wind_speed_image.lt(wind_speed_threshold)
                
                # Update mask to filter the image
                return image.updateMask(mask)

            masked_collection = image_collection_with_wind.map(wind_mask)

            def count_valid_pixels(image):
                non_zero_pixels = image.select(0).reduceRegion(
                    reducer=ee.Reducer.count(),
                    geometry=roi,
                    scale=30,
                    maxPixels=1e9
                ).get('VV')  # Remplacer 'VV' par le nom d'une bande qui existe dans votre image (VV ou VH)

                return image.set('non_zero_pixels', non_zero_pixels)

            counted_collection = masked_collection.map(count_valid_pixels)

            def filter_by_pixel_count(image):
                return ee.Image(image).set('mask', ee.Number(image.get('non_zero_pixels')).gt(min_pixels))

            filtered_collection = counted_collection.map(filter_by_pixel_count)
            final_collection = filtered_collection.filterMetadata('mask', 'equals', 1).map(lambda image : image.clip(roi))

            return final_collection
        
        # Filtrage des images pour la collection descendante
        wind_filtered_vv_vh_iw_desc = filter_images_by_wind_speed(vv_vh_iw_desc, roi) 

        # Filtrage des images pour la collection ascendante
        wind_filtered_vv_vh_iw_asc = filter_images_by_wind_speed(vv_vh_iw_asc,roi) 

    #### NO WIND FILTER #################################################################################
    else : 
        
        def count_valid_pixels(image):
                non_zero_pixels = image.select(0).reduceRegion(
                    reducer=ee.Reducer.count(),
                    geometry=roi,
                    scale=30,
                    maxPixels=1e9
                ).get('VV')  # Remplacer 'VV' par le nom d'une bande qui existe dans votre image (VV ou VH)

                return image.set('non_zero_pixels', non_zero_pixels)

        pixel_counted_vv_vh_iw_desc = vv_vh_iw_desc.map(count_valid_pixels)
        pixel_counted_vv_vh_iw_asc = vv_vh_iw_asc.map(count_valid_pixels)

        def filter_by_pixel_count(image,min_pixels=1000):
            return ee.Image(image).set('mask', ee.Number(image.get('non_zero_pixels')).gt(min_pixels))

        add_properties_vv_vh_iw_desc = pixel_counted_vv_vh_iw_desc.map(filter_by_pixel_count)
        wind_filtered_vv_vh_iw_desc = add_properties_vv_vh_iw_desc.filterMetadata('mask', 'equals', 1).map(lambda image : image.clip(roi))
        add_properties_vv_vh_iw_asc = pixel_counted_vv_vh_iw_asc.map(filter_by_pixel_count)
        wind_filtered_vv_vh_iw_asc = add_properties_vv_vh_iw_asc.filterMetadata('mask', 'equals', 1).map(lambda image : image.clip(roi))
    
    #### ANGLE CORRECTION #################################################################################

    def toGammaVV_VH(image):
        # Correction for VV band
        corr_VV = image.select('VV').subtract(
            image.select('angle').multiply(math.pi/180.0).cos().log10().multiply(10.0)
        ).rename('VV')

        # Correction for VH band
        corr_VH = image.select('VH').subtract(
            image.select('angle').multiply(math.pi/180.0).cos().log10().multiply(10.0)
        ).rename('VH')

        # Add the corrected bands and the angle band to a new image
        return image.addBands([corr_VV, corr_VH], overwrite=True)

    angle_corrected_vv_vh_iw_desc = wind_filtered_vv_vh_iw_desc.map(toGammaVV_VH)
    angle_corrected_vv_vh_iw_asc = wind_filtered_vv_vh_iw_asc.map(toGammaVV_VH)

    #### COMPOSITE CREATION #################################################################################
    
    def composite_12day_intervals(collection, start_date, end_date):
    
        def wrap(date):
            date = ee.Date(date)
            interval_collection = collection.filterDate(date, date.advance(12, 'day'))
            return interval_collection.median().set('system:time_start', date.millis())
        
        days = ee.List.sequence(0, end_date.difference(start_date, 'day').subtract(1), 12)
        dates = days.map(lambda day: start_date.advance(day, 'day'))
        
        col = ee.ImageCollection(dates.map(wrap))
        
        def format_date(image):
            return image.set('date', image.date().format('Y-MM-dd'))
        
        col = col.map(format_date)
        
        return col
    
    def format_date(image):
            return image.set('date', image.date().format('Y-MM-dd'))
    
    def filter_empty_bands(image):
        band_count = image.bandNames().size()
        return image.set('band_count', band_count)

    monthly_vv_vh_iw_asc = angle_corrected_vv_vh_iw_asc.map(format_date)
    monthly_vv_vh_iw_desc = angle_corrected_vv_vh_iw_desc.map(format_date)
    #monthly_vv_vh_iw_asc = composite_12day_intervals(angle_corrected_vv_vh_iw_asc, start_date, end_date)
    #monthly_vv_vh_iw_desc = composite_12day_intervals(angle_corrected_vv_vh_iw_desc, start_date, end_date)
    monthly_vv_vh_iw_asc = monthly_vv_vh_iw_asc.map(filter_empty_bands).filter(ee.Filter.eq('band_count', 3))
    monthly_vv_vh_iw_desc = monthly_vv_vh_iw_desc.map(filter_empty_bands).filter(ee.Filter.eq('band_count', 3))

    #### LEE FILTER #################################################################################

    def toNatural(img):
        """Function to convert from dB"""
        return ee.Image(10.0).pow(img.select(0).divide(10.0))

    def toDB(img):
        """Function to convert to dB"""
        return ee.Image(img).log10().multiply(10.0)

    def RefinedLee(img):
        """The RL speckle filter
        img must be in natural units, i.e. not in dB!
        Set up 3x3 kernels"""
        bandNames = img.bandNames()
        img = toNatural(img)
        
        weights3 = ee.List.repeat(ee.List.repeat(1,3),3)
        kernel3 = ee.Kernel.fixed(3,3, weights3, 1, 1, False)

        mean3 = img.reduceNeighborhood(ee.Reducer.mean(), kernel3)
        variance3 = img.reduceNeighborhood(ee.Reducer.variance(), kernel3)

        # Use a sample of the 3x3 windows inside a 7x7 windows to determine gradients and directions
        sample_weights = ee.List([[0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0], [0,1,0,1,0,1,0], [0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0]])

        sample_kernel = ee.Kernel.fixed(7,7, sample_weights, 3,3, False)

        # Calculate mean and variance for the sampled windows and store as 9 bands
        sample_mean = mean3.neighborhoodToBands(sample_kernel)
        sample_var = variance3.neighborhoodToBands(sample_kernel)

        # Determine the 4 gradients for the sampled windows
        gradients = sample_mean.select(1).subtract(sample_mean.select(7)).abs()
        gradients = gradients.addBands(sample_mean.select(6).subtract(sample_mean.select(2)).abs())
        gradients = gradients.addBands(sample_mean.select(3).subtract(sample_mean.select(5)).abs())
        gradients = gradients.addBands(sample_mean.select(0).subtract(sample_mean.select(8)).abs())

        # And find the maximum gradient amongst gradient bands
        max_gradient = gradients.reduce(ee.Reducer.max())

        # Create a mask for band pixels that are the maximum gradient
        gradmask = gradients.eq(max_gradient)

        # duplicate gradmask bands: each gradient represents 2 directions
        gradmask = gradmask.addBands(gradmask)

        # Determine the 8 directions
        directions = sample_mean.select(1).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(7))).multiply(1)
        directions = directions.addBands(sample_mean.select(6).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(2))).multiply(2))
        directions = directions.addBands(sample_mean.select(3).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(5))).multiply(3))
        directions = directions.addBands(sample_mean.select(0).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(8))).multiply(4))
        # The next 4 are the not() of the previous 4
        directions = directions.addBands(directions.select(0).Not().multiply(5))
        directions = directions.addBands(directions.select(1).Not().multiply(6))
        directions = directions.addBands(directions.select(2).Not().multiply(7))
        directions = directions.addBands(directions.select(3).Not().multiply(8))

        # Mask all values that are not 1-8
        directions = directions.updateMask(gradmask)

        # "collapse" the stack into a singe band image (due to masking, each pixel has just one value (1-8) in it's directional band, and is otherwise masked)
        directions = directions.reduce(ee.Reducer.sum())

        #pal = ['ffffff','ff0000','ffff00', '00ff00', '00ffff', '0000ff', 'ff00ff', '000000']
        #Map.addLayer(directions.reduce(ee.Reducer.sum()), {min:1, max:8, palette: pal}, 'Directions', False)

        sample_stats = sample_var.divide(sample_mean.multiply(sample_mean))

        # Calculate localNoiseVariance
        sigmaV = sample_stats.toArray().arraySort().arraySlice(0,0,5).arrayReduce(ee.Reducer.mean(), [0])

        # Set up the 7*7 kernels for directional statistics
        rect_weights = ee.List.repeat(ee.List.repeat(0,7),3).cat(ee.List.repeat(ee.List.repeat(1,7),4))

        diag_weights = ee.List([[1,0,0,0,0,0,0], [1,1,0,0,0,0,0], [1,1,1,0,0,0,0],
        [1,1,1,1,0,0,0], [1,1,1,1,1,0,0], [1,1,1,1,1,1,0], [1,1,1,1,1,1,1]])

        rect_kernel = ee.Kernel.fixed(7,7, rect_weights, 3, 3, False)
        diag_kernel = ee.Kernel.fixed(7,7, diag_weights, 3, 3, False)

        # Create stacks for mean and variance using the original kernels. Mask with relevant direction.
        dir_mean = img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel).updateMask(directions.eq(1))
        dir_var = img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel).updateMask(directions.eq(1))

        dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel).updateMask(directions.eq(2)))
        dir_var= dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel).updateMask(directions.eq(2)))

        # and add the bands for rotated kernels
        for i in range(1,4):
            dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)))
            dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)))
            dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)))
            dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)))
            
        # "collapse" the stack into a single band image (due to masking, each pixel has just one value in it's directional band, and is otherwise masked)
        dir_mean = dir_mean.reduce(ee.Reducer.sum())
        dir_var = dir_var.reduce(ee.Reducer.sum())

        # A finally generate the filtered value
        varX = dir_var.subtract(dir_mean.multiply(dir_mean).multiply(sigmaV)).divide(sigmaV.add(1.0))

        b = varX.divide(dir_var)

        result = dir_mean.add(b.multiply(img.subtract(dir_mean))).arrayFlatten([['sum']]) \
                .float()
                
        return ee.Image(toDB(result).rename('leeFilter')).set('date',img.get('date'))#.rename(bandNames)

    lee_filtered_vv_vh_iw_desc = monthly_vv_vh_iw_desc.map(RefinedLee)
    lee_filtered_vv_vh_iw_asc = monthly_vv_vh_iw_asc.map(RefinedLee)

    def addDateProperty(image, originalCollection):
        index = image.get('system:index')
        original_image = ee.Image(originalCollection.filterMetadata('system:index', 'equals', index).first())
        date = original_image.get('date')
        return image.set('date', date)

    # Ajoute la propriété de date à la collection descendante
    lee_filtered_vv_vh_iw_desc = lee_filtered_vv_vh_iw_desc.map(
    lambda img: addDateProperty(img, monthly_vv_vh_iw_desc)
    )

    # Ajoute la propriété de date à la collection ascendante
    lee_filtered_vv_vh_iw_asc = lee_filtered_vv_vh_iw_asc.map(
    lambda img: addDateProperty(img, monthly_vv_vh_iw_asc)
    )
    #### CLASSIFIER #################################################################################

    def radarClassifier (input, region, scale, clusterNumber):
    # Make the training dataset.
        training = input.sample(region = region,scale = scale,numPixels = 10000) #geometry ??
        # Instantiate the clusterer and train it.
        # Weka kmeans clusterer
        #clusterer = ee.Clusterer.wekaKMeans(clusterNumber, 1).train(training)
        clusterer = ee.Clusterer.wekaKMeans(clusterNumber, distanceFunction="Euclidean", maxIterations=75).train(training)
        # Cluster the input using the trained clusterer.
        # Display the clusters with random colors.
        return input.addBands(input.cluster(clusterer)) #|| undefined

    # Filtrage des images pour la collection descendante
    clustered_vv_vh_iw_desc = lee_filtered_vv_vh_iw_desc.map(lambda image : radarClassifier(image, roi, 10, 2))
    # Filtrage des images pour la collection ascendante
    clustered_vv_vh_iw_asc = lee_filtered_vv_vh_iw_asc.map(lambda image : radarClassifier(image, roi, 10, 2))

    #### WATER IDENTIFICATION #################################################################################

    def identify_water_cluster(clustered_image):
        # Calculer la moyenne des bandes spécifiées pour chaque cluster
        cluster_0 = clustered_image.updateMask(clustered_image.select('cluster').eq(0))
        cluster_1 = clustered_image.updateMask(clustered_image.select('cluster').eq(1))

        # Réduire la région pour calculer la moyenne des bandes VV et VH pour chaque cluster
        cluster_0_mean = cluster_0.select('leeFilter').reduceRegion(reducer=ee.Reducer.mean(), geometry=roi, scale=10)
        cluster_1_mean = cluster_1.select('leeFilter').reduceRegion(reducer=ee.Reducer.mean(), geometry=roi, scale=10)

        # Extraire les moyennes de VV et VH pour chaque cluster
        mean_vv_cluster_0 = cluster_0_mean.getNumber('leeFilter')
        mean_vv_cluster_1 = cluster_1_mean.getNumber('leeFilter')

        # Identifier le cluster avec la somme la plus petite de VV et VH
        water_cluster = ee.Image(ee.Algorithms.If(
            mean_vv_cluster_0.lt(mean_vv_cluster_1),
            clustered_image.select('cluster').eq(0),
            clustered_image.select('cluster').eq(1)
        )).rename('water')

        # Créer une image binaire avec 1 pour l'eau et 0 sinon
        water_image = water_cluster.unmask(0).rename('water').set('date', clustered_image.get('date'))

        return water_image

    # Appliquer la fonction à vos collections d'images
    water_cluster_vv_vh_iw_desc = clustered_vv_vh_iw_desc.map(identify_water_cluster)
    water_cluster_vv_vh_iw_asc = clustered_vv_vh_iw_asc.map(identify_water_cluster)

    return water_cluster_vv_vh_iw_desc, water_cluster_vv_vh_iw_asc
