# -*- coding: utf-8 -*-
"""
Author: Ty Nietupski (ty.nietupski@oregonstate.edu)

The MIT License

Copyright © 2021 Ty Nietupski

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""

import ee


def calcSpecDist(maskedLandsat, modSorted_t01, numPixels, pixPositions):
    """
    Calculate the mean absolute spectral distance between the Landsat and\
    MODIS pixels.

    Parameters
    ----------
    maskedLandsat : ee_list.List
        Landsat neighborhood images masked and prepared with the prepLandsat
        function.
    modSorted_t01 : ee_list.List
        MODIS neighborhood images prepared with the prepMODIS function.
    numPixels : ee_number.Number
        Total number of pixels in the kernel.
    pixPositions : ee_list.List
        Names for each index in the window (e.g., "_0_0", "_0_1" ...).

    Returns
    -------
    image.Image
        Image with the spectral distance between each pixel in neighborhood.

    """
    sDist = ee.List.sequence(0, numPixels.subtract(1)) \
        .map(lambda index:
             ee.Image(maskedLandsat.get(index)) \
             .subtract(ee.Image(modSorted_t01.get(index))) \
             .abs() \
             .reduce(ee.Reducer.mean()))

    return ee.ImageCollection(sDist) \
        .toBands() \
        .rename(pixPositions.map(lambda name:
                                 ee.String('sDist').cat(name)))


def calcSpatDist(positions):
    """
    Calculate the spatial distance between each pixel in the window and the\
    central pixel.

    Parameters
    ----------
    positions : ee_list.List
        Names for each index in the window (e.g., "_0_0", "_0_1" ...).

    Returns
    -------
    image.Image
        Image with the spatial distance between each pixel in neighborhood.

    """
    # window width
    w2 = positions.length().sqrt().subtract(1).divide(2)

    # distance to each pixel in window
    dist = positions.map(lambda position:
                         ee.Image.constant(ee.Number(1)
                                           .add(ee.Number.parse(
                                               ee.String(position) \
                                               .match('(-?[0-9]+)', 'g')\
                                               .get(0)) \
                                           .pow(2) \
                                           .add(ee.Number.parse(
                                               ee.String(position) \
                                               .match('(-?[0-9]+)', 'g')\
                                               .get(1))\
                                           .pow(2)) \
                                           .sqrt() \
                                           .divide(w2))))

    return ee.ImageCollection(dist) \
        .toBands() \
        .rename(positions.map(lambda bn: ee.String('corr').cat(bn)))

def calcWeight(spatDist, specDist):
    """
    Create diagonal weight matrix from a combination of the spatial and\
    spectral distance images.

    Parameters
    ----------
    spatDist : image.Image
        Spatial distance to each pixel in the kernel (window).
    specDist : image.Image
        Spectral distance to each pixel in the kernel (window).

    Returns
    -------
    image.Image (array image)
        Weight for all similar pixels in the window.

    """
    disIndex = specDist.multiply(spatDist)
    num = ee.Image.constant(1).divide(disIndex)
    sumNum = ee.ImageCollection(num.bandNames()
                                .map(lambda bn:
                                     num.select([bn]).rename('num'))) \
        .sum()
    W = num.divide(sumNum)

    return W.unmask().toArray().toArray(1).matrixToDiag()


def calcConversionCoeff(maskedLandsat, modSorted_t01,
                        doys, numPixels, commonBandNames):
    """
    Perform linear regression between all Landsat and MODIS pixels in the\
    window.

    Parameters
    ----------
    maskedLandsat : ee_list.List
        Landsat neighborhood images masked and prepared with the prepLandsat
        function.
    modSorted_t01 : ee_list.List
        MODIS neighborhood images prepared with the prepMODIS function.
    doys : ee_list.List
        List of day of year associated with t0 and t1.
    numPixels : ee_number.Number
        Total number of pixels in the kernel.
    commonBandNames : ee_list.List
        Names of bands to use in fusion.

    Returns
    -------
    coeffs : image.Image (array image)
        Scaling coefficients for the window.

    """
    # reformat the landsat & modis data
    lanMod = doys \
        .map(lambda doy:
             ee.List.sequence(0, numPixels.subtract(1)) \
             .map(lambda index:
                  ee.Image.constant(1).rename(['intercept']) \
                  .addBands(ee.Image(modSorted_t01.get(index))\
                                .select(ee.String(doy).cat('.+')) \
                            .rename(commonBandNames\
                                    .map(lambda bn:
                                         ee.String(bn).cat('_modis')))) \
                  .addBands(ee.Image(maskedLandsat.get(index))\
                                .select(ee.String(doy).cat('.+')) \
                            .rename(commonBandNames\
                                    .map(lambda bn:
                                         ee.String(bn).cat('_landsat'))))))

    # when we convert this collection to an array we get a 2-D array image
    # where the 0 axis is the similar pixel images and the 1 axis is the bands
    # of landsat and modis
    lanMod = ee.ImageCollection(lanMod.flatten())

    # solve for conversion coefficients using linear regression reducer
    coeffs = lanMod \
        .reduce(ee.Reducer.linearRegression(commonBandNames.length().add(1),
                                            commonBandNames.length()))\
        .select([0], ['coefficients']) \
        .arraySlice(0, 1, commonBandNames.length().add(1))

    return coeffs


def predictLandsat(landsat_t01, modSorted_t01,
                   doys, modSorted_tp, weights,
                   coeffs, commonBandNames, numPixels):
    """
    Use weights and coefficients to predict landsat from MODIS. We predict\
    from t0 and t1 and then merge predictions, weighting by the temporal\
    proximity of predicted image to t0 and t1, to get the final prediction.

    Parameters
    ----------
    landsat_t01 : ee_list.List
        Landsat images at time 0 (t0) and time 1 (t1).
    modSorted_t01 : ee_list.List
        MODIS neighborhood images prepared with the prepMODIS function.
    doys : ee_list.List
        List of day of year associated with t0 and t1.
    modSorted_tp : ee_list.List
        MODIS neighborhood image for the prediction date prepared with
        prepMODIS function.
    weights : image.Image (array image)
        Weight for all similar pixels in the window.
    coeffs : image.Image (array image)
        Scaling coefficients for the window.
    commonBandNames : ee_list.List
        Names of bands to use in fusion.
    numPixels : ee_number.Number
        Total number of pixels in the kernel.

    Returns
    -------
    image.Image
        Predicted Landsat image for date corresponding to modSorted_tp.

    """
    # split apart the modis images based on the doy
 
