# This script runs the legacy SVM against all five data sets,
# computing a 10 fold cross validation error computation and printing
# the data out to the screen.
import numpy as numpy

from sklearn import svm
from sklearn.model_selection import cross_val_score

# Defines the readData function
execfile('readData.py')
    
fileNames = ['QuietData.csv',
             'UnalignedData.csv',
             'BoundaryData.csv',
             'BeamPatternData.csv',
             'BeamShallowData.csv']

for index in range(len(fileNames)):
    fileName = fileNames[index]

    # First read the quiet data and do a quick validation on it.
    parsedData = readData(fileName)
    model = svm.SVC(kernel='rbf',gamma=0.001, C=100)
    print fileName," SVM(0.001,100) 10 Fold CVE: ",
    score = cross_val_score(model, parsedData.x, y=parsedData.y, cv=10);
    print numpy.min(score)," < ",numpy.mean(score)," < ",numpy.max(score)

