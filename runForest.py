# This script runs the random forest classifier using default parameters
# on all of the data sets, computing a 10 fold CVE and printing it to the
# screen
import numpy as numpy
from sklearn.model_selection import cross_val_score
from sklearn.ensemble import RandomForestClassifier


# Defines the readData function
execfile('readData.py')

# Define all the data sets
fileNames = ['QuietData.csv',
             'UnalignedData.csv',
             'BoundaryData.csv',
             'BeamPatternData.csv',
             'BeamShallowData.csv']

for index in range(len(fileNames)):
    parsedData = readData(fileNames[index])
    model = RandomForestClassifier()
    print fileNames[index],": 10 Fold CVE RandomForest ",
    score = cross_val_score(model, parsedData.x, y=parsedData.y,cv=10);
    print numpy.min(score)," < ",numpy.mean(score)," < ",numpy.max(score)

    
