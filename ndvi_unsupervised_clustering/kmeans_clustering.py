# -*- coding: utf-8 -*- 

""" 
Author:  Magnier Morgane (morgane.magnier@vattenfall.com) 

Copyright © 2024 Magnier Morgane
This work is part of a thesis project. The copyright of the thesis itself belongs to the student Magnier Morgane. 
Vattenfall has the right to use the thesis work, including its findings, methods, and conclusions, in its operations. Any material generated within the framework of this thesis that is subject to intellectual property protection (e.g., source code, computer program, design, or invention) is the property of Vattenfall, unless otherwise agreed in writing. Permission is hereby granted to view, copy, and share this work for educational or personal purposes only, provided that this copyright notice and this permission notice appear in all copies.
"""

import ee, geemap

####################################### Kmeans clustering #######################################

def kmeans_clustering_S2(image, roi, bands = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B9', 'B11', 'B12'], scale = 10, numPixels = 1500, nClusters = 4,init = 2, distanceFunction = "Euclidean", maxIterations = 75): 
    
    training = image.select(bands).sample(
    region = roi,
    scale = scale,
    numPixels = numPixels,
    geometries=True
    )

    kmeans = ee.Clusterer.wekaKMeans(nClusters = nClusters, init = init, distanceFunction = distanceFunction, maxIterations = maxIterations).train(training)
    kmean_rs = image.cluster(kmeans)

    return kmean_rs

####################################### CLUSTER IDENTIFICATION IF nClusters == 4 #######################################

def instentiate_collection(nClusters):

    cluster_numbers = ee.List.sequence(0, nClusters - 1)
    cluster_images = cluster_numbers.map(lambda num_cluster: ee.Image.constant(num_cluster).rename('cluster').set({'num_cluster': num_cluster}))
    cluster_images_collection = ee.ImageCollection(cluster_images)

    return cluster_images_collection

def mask_image_by_cluster(img, clusters, num_cluster):

    cluster_mask = clusters.eq(num_cluster)
    masked_img = img.updateMask(cluster_mask)
    return masked_img.set({'num_cluster' : num_cluster.get('num_cluster')})

def calculate_median(img, roi):

    median_values = img.reduceRegion(reducer=ee.Reducer.median(), geometry = roi, scale=10)
    
    feature = ee.Feature(None, {
        'median_NDWI': median_values.get('NDWI'),
        'median_NDVI': median_values.get('NDVI'),
        'median_NDGI': median_values.get('NDGI'),
        'median_BI': median_values.get('BI'),
        'num_cluster': img.get('num_cluster')
    })
    
    return feature

def identify_clusters(img, roi, clusters, nClusters = 4): 

    cluster_images_collection = instentiate_collection(nClusters)
    masked_images_collection = cluster_images_collection.map(lambda num_cluster : mask_image_by_cluster(img, clusters, num_cluster))
    cluster_median_features = masked_images_collection.map(lambda image : calculate_median(image, roi))

    cluster_median_features_df = ee.data.computeFeatures({
        'expression': cluster_median_features,
        'fileFormat': 'PANDAS_DATAFRAME'
        })

    cluster_median_features_df_2 = cluster_median_features_df.copy()

    max_ndgi_cluster = int(cluster_median_features_df_2.loc[cluster_median_features_df_2['median_NDGI'].idxmax(), 'num_cluster'])
    cluster_median_features_df_2 = cluster_median_features_df_2[cluster_median_features_df_2['num_cluster'] != max_ndgi_cluster]
    max_ndvi_cluster = int(cluster_median_features_df_2.loc[cluster_median_features_df_2['median_NDVI'].idxmax(), 'num_cluster'])
    max_ndwi_cluster = int(cluster_median_features_df_2.loc[cluster_median_features_df_2['median_NDWI'].idxmax(), 'num_cluster'])
    cluster_median_features_df_2 = cluster_median_features_df_2[cluster_median_features_df_2['num_cluster'] != max_ndwi_cluster]
    max_bi_cluster = int(cluster_median_features_df_2.loc[cluster_median_features_df_2['median_BI'].idxmax(), 'num_cluster'])
    
    clusters_annotated = (clusters.where(clusters.eq(max_ndwi_cluster), 0).where(clusters.eq(max_ndvi_cluster), 1).where(clusters.eq(max_ndgi_cluster), 2).where(clusters.eq(max_bi_cluster), 3))

    return clusters_annotated

def cut_img_by_clusters(img, clusters_bands): 

    water = img.updateMask(clusters_bands.select('0_water'))
    grass = img.updateMask(clusters_bands.select('1_trees'))
    trees = img.updateMask(clusters_bands.select('2_grass'))
    bare_land = img.updateMask(clusters_bands.select('3_bare'))

    return water, grass, trees, bare_land 

####################################### Visualization #######################################

def plot_clusters_in_wetlands_S2(img, roi, clusters_annotated, wetlands):

    rgbVisS2 = {'bands': ['B4',  'B3',  'B2'], 'min': 0, 'max':  0.1}

    waterVis = {'min': 0, 'max': 1, 'palette': ['#0000FF']}
    bareVis = {'min': 0, 'max': 1, 'palette': ['#FFFF00']}
    treesVis = {'min': 0, 'max': 1, 'palette': ['#008000']}
    grassVis = {'min': 0, 'max': 1, 'palette': ['#800080']}

    m = geemap.Map()
    m.centerObject(roi, 14)

    m.addLayer(img.scale().clip(roi),rgbVisS2,'img')
    m.addLayer(clusters_annotated.eq(0).selfMask(), waterVis, 'Water')
    m.addLayer(clusters_annotated.eq(1).selfMask(), treesVis, 'Trees')
    m.addLayer(clusters_annotated.eq(2).selfMask(), grassVis, 'Grass')
    m.addLayer(clusters_annotated.eq(3).selfMask(), bareVis, 'Bare')
    wetlands_contours = wetlands.subtract(wetlands.focalMin(20, 'square', 'meters'))
    m.addLayer(wetlands_contours.updateMask(wetlands_contours),{'min': 0, 'max': 1, 'palette': ['blue']})

    labels = ['Water', 'Bare/Small grass', 'Trees', 'Grass']
    colors = ['#0000FF', '#FFFF00', '#008000', '#800080']

    m.add_legend(labels=labels, colors=colors, position='bottomright')

    return m

def plot_clusters_S2(img, roi, clusters_bands, wetlands):

    rgbVisS2 = {'bands': ['B4',  'B3',  'B2'], 'min': 0, 'max':  0.1}

    waterVis = {'min': 0, 'max': 3, 'palette': ['#0000FF']}
    bareVis = {'min': 0, 'max': 3, 'palette': ['#FFFF00']}
    treesVis = {'min': 0, 'max': 3, 'palette': ['#008000']}
    grassVis = {'min': 0, 'max': 3, 'palette': ['#800080']}

    m = geemap.Map()
    m.centerObject(roi, 14)

    m.addLayer(img.scale().clip(roi),rgbVisS2,'img')
    m.addLayer(clusters_bands.select('0_water').clip(roi), waterVis, 'Water')
    m.addLayer(clusters_bands.select('1_trees').clip(roi), treesVis, 'Trees')
    m.addLayer(clusters_bands.select('2_grass').clip(roi), grassVis, 'Grass')
    m.addLayer(clusters_bands.select('3_bare').clip(roi), bareVis, 'Bare')
    wetlands_contours = wetlands.subtract(wetlands.focalMin(20, 'square', 'meters'))
    m.addLayer(wetlands_contours.updateMask(wetlands_contours),{'min': 0, 'max': 1, 'palette': ['blue']})
    labels = ['Water', 'Bare', 'Trees', 'Grass']
    colors = ['#0000FF', '#FFFF00', '#008000', '#800080']

    m.add_legend(labels=labels, colors=colors, position='bottomright')

    return m


