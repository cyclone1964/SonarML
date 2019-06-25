# This script runs the models using different amounts of pre- and post-
# frames on the unaligned data set
import numpy as numpy
from sklearn import svm
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.neural_network import MLPClassifier

# Define the function for reading data
execfile('readData.py')
    
# Define the names of the datsets
fileNames = ['QuietData.csv',
             'UnalignedData.csv',
             'BoundaryData.csv',
             'BeamPatternData.csv',
             'BeamShallowData.csv']

# Run for each radius ..
for radius in range(5):
    # .. and for each file.
    for fileIndex in range(len(fileNames)):
        # Read the file
        parsedData = readData(fileNames[fileIndex])

        # If the radius is non-zero we need to form the agglomerated
        # stuff, which is a big mess to get rid of the frames that
        # don't have support. As currently coded, this includes the
        # annotations.
        if (radius > 0):

            # Get the labels
            y = parsedData.y[radius:-radius]
            x = parsedData.annotations[radius:-radius,:]

            # do the first ones, but we have to do the last one
            # special because in python, the derferenceing [a:0] is
            # different than [a:]
            for index in range(2*radius):
                x = numpy.concatenate((x,
                                       parsedData.x[index:-(2*radius-index),:]),
                                      axis=1);
            x = numpy.concatenate((x,
                                   parsedData.x[2*radius:]),
                                  axis=1);
            # Now form a list of the bad frames, which are every 32
            # starting from the ones in the radius.
            invalid = []
            for i in range(0,parsedData.numCycles-1):
                for j in range(31-radius,31+radius):
                    invalid.append(j + 32 * i)
            invalid.sort()

            # Now for the valid ones with the worst piece of code ever
            valid = []
            for index in range(len(y)):
                if (not (index in invalid)):
                    valid.append(index)

            y = y[valid]
            x = x[valid,:];
        else:
            # If ther is no radius, just use the data as it sits. 
            y = parsedData.y
            x = numpy.concatenate((parsedData.annotations,
                                       parsedData.x),
                                       axis=1);

        # Now lets normalize the input data
        dim = numpy.shape(x)

        for index in range(dim[0]):
            limits = [numpy.min(x[index,4:]),
                      numpy.max(x[index,4:])]
            x[index,4:] = (x[index,4:] - limits[0])/(limits[1]-limits[0])

        model = svm.SVC(kernel='rbf',gamma=0.001, C=100)
        print fileNames[fileIndex],": 10 Fold CVE SVM ",radius," ",
        score = cross_val_score(model, x, y=y,cv=10);
        print numpy.min(score)," < ",
        print numpy.mean(score)," < ",numpy.max(score)
            
        model = RandomForestClassifier()
        print fileNames[fileIndex],": 10 Fold CVE RF ",radius," ",
        score = cross_val_score(model, x, y=y,cv=10);
        print numpy.min(score)," < ",
        print numpy.mean(score)," < ",numpy.max(score)

        model = MLPClassifier(activation='logistic',
                              solver='adam',
                              hidden_layer_sizes = (64))
            
        print fileNames[fileIndex],": 10 Fold CVE ANN ",radius," ",
        score = cross_val_score(model, x, y=y,cv=10);
        print numpy.min(score)," < ",
        print numpy.mean(score)," < ",numpy.max(score)

            
                                  
