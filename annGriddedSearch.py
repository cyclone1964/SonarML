# This script does a gridded search on neural network architectures
# to find the best one for the legacy Quiet data set
import numpy as numpy
from sklearn.model_selection import cross_val_score
from sklearn.neural_network import MLPClassifier

# Defines the readData function
execfile('readData.py')

fileName = 'QuietData.csv'
layerSets = [(4), (8), (16), (32), (64), (4,4),
             (8,8), (16,16), (32,32), (64,64)]

# do each layer set, which is to say architectures
for setIndex in range(len(layerSets)):

    layerSet = layerSets[setIndex]
    parsedData = readData(fileName)
    model = MLPClassifier(activation='logistic',
                          solver='lbfgs',
                          hidden_layer_sizes = layerSet)
    print fileName," ANN(sigmoid,",layerSet,") 10 Fold CVE: ",
    score = cross_val_score(model, parsedData.x, y=parsedData.y,cv=10);
    print numpy.min(score)," < ",numpy.mean(score)," < ",numpy.max(score)

