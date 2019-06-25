# This script runs the neural net using the optimal architecture found
# by the gridded search on all of the data sets, computing a 10 fold
# CVE and printing it to the screen
import numpy as numpy

from sklearn.model_selection import cross_val_score
from sklearn.neural_network import MLPClassifier


# Defines the readData function
execfile('readData.py')

# First read the quiet data and do a quick validation on it.
fileNames = ['QuietData.csv',
             'UnalignedData.csv',
             'BoundaryData.csv',
             'BeamPatternData.csv',
             'BeamShallowData.csv']
layerSet = (64)

for index in range(len(fileNames)):
    parsedData = readData(fileNames[index])
    model = MLPClassifier(activation='logistic',
                          solver='lbfgs',
                          hidden_layer_sizes = layerSet)
    print fileNames[index],": 10 Fold CVE RandomForest ",
    score = cross_val_score(model, parsedData.x, y=parsedData.y,cv=10);
    print numpy.min(score)," < ",numpy.mean(score)," < ",numpy.max(score)

    
