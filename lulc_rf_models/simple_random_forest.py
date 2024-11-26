

# -*- coding: utf-8 -*-

""" 
Author:  Magnier Morgane (morgane.magnier@vattenfall.com) 

Copyright © 2024 Magnier Morgane
This work is part of a thesis project. The copyright of the thesis itself belongs to the student Magnier Morgane. 
Vattenfall has the right to use the thesis work, including its findings, methods, and conclusions, in its operations. Any material generated within the framework of this thesis that is subject to intellectual property protection (e.g., source code, computer program, design, or invention) is the property of Vattenfall, unless otherwise agreed in writing. Permission is hereby granted to
view, copy, and share this work for educational or personal purposes only, provided that this copyright notice and this permission notice appear in all copies. 
"""


import ee 

def Random_Forest(composite, polygons, numberOfTrees = 500, variablesPerSplit = 2, minLeafPopulation = 1, bagFraction = 0.5, maxNodes = 7, seed = 11): 

    training_data = composite.sampleRegions(**{
    'collection': polygons,
    'properties': ['class'],
    'scale': 10
    })

    # Training & validation splitting
    sample = training_data.randomColumn()
    train_sample = sample.filter('random <= 0.8')
    val_sample = sample.filter('random > 0.8')

    clsfr_rf = ee.Classifier.smileRandomForest(numberOfTrees, variablesPerSplit, minLeafPopulation ,, bagFraction, maxNodes, seed)

    trained_clsfr_rf = clsfr_rf.train(
        features=train_sample,
        classProperty='class',
        inputProperties=composite.bandNames(),)

    display('Results of trained classifier RF', trained_clsfr_rf.explain())

    train_accuracy_rf = trained_clsfr_rf.confusionMatrix()

    display('Training error matrix', train_accuracy_rf) #No display (client side op) for mapping in a collection
    display('Training overall accuracy', train_accuracy_rf.accuracy())

    val_sample_rf = val_sample.classify(trained_clsfr_rf)

    val_accuracy_rf = val_sample_rf.errorMatrix('class', 'classification')

    display('Validation error matrix', val_accuracy_rf)
    display('Validation accuracy', val_accuracy_rf.accuracy())

    return trained_clsfr_rf
