# This is an attempt to apply an SVM to the data. Specficially, we
# train 4 SVMS each 48 bins wide plus annotations.
import sys
import string
import matplotlib.pyplot as pyplot
import numpy as numpy
from sklearn import svm
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score

# This holds the data structure
class ParsedData:
    def __init__(self):
        self.x = []
        self.y = []
        self.cycles = []
        self.frames = []
        self.dopplers = []
        self.numRows = 0
        self.numColums = 0

# This function reads in the data file. 
def readData(fileName):

    # We read in the file. This is a CSV with the following columns
    #
    # Target/NoTarget
    # TargetDoppler,
    # CycleIndex,
    # FrameIndex,
    # FrameRange,
    # PlatformSpeed
    # Azimuth of Beam
    # Elevation of Beam
    # 128 Bins
    print "Read File <",fileName,">"
    fileID = open(fileName,'r')
    lines = fileID.readlines()
    fileID.close();

    # Compute the number of rows and columns
    output = ParsedData()
    output.numRows = len(lines)
    fields = lines[1].split(",")
    output.numColumns = len(fields)-4
    print "    Reading ",output.numRows," rows of ", output.numColumns, "fields"

    # Initialize the label and the data array
    output.x = numpy.zeros([output.numRows, output.numColumns])
    output.y = numpy.zeros([output.numRows])
    output.frames = numpy.zeros([output.numRows])
    output.cycles = numpy.zeros([output.numRows])
    output.dopplers = numpy.zeros([output.numRows])

    for lineIndex in range(output.numRows):
        line = lines[lineIndex]
        fields = line.split(",")
        output.y[lineIndex] = fields[0]
        output.dopplers[lineIndex] = fields[1]
        output.cycles[lineIndex] = fields[2]
        output.frames[lineIndex] = fields[3]
        output.x[lineIndex,:] = fields[4:]
    print "    Done!"
    
    # Now, for fun, let's try normalizing these
    #limits = [numpy.min(output.x[:,8:]), numpy.max(output.x[:,8:])]
    #output.x[:,8:] = (output.x[:,8:] - limits[0])/(limits[1] - limits[0])
    
    
    return output

# The main function. Right now the steps of things I intend to do for
# this project and which this function does is:
#
# Re-run the legacy problem using the new simulation and re-producing
# the results.
#
# Process the advanced data set with reverberation of all types and
# check the performance using an SVM with the same parameters
#
# 
#
# 
if __name__ == "__main__":

    # First read the quiet data and do a quick validation on it.
    parsedData = readData('QuietData.csv')
    model = svm.SVC(kernel='rbf',gamma=0.001, C=100)
    print "    10 Fold CVE C=100, Gamma = 0.01 ..."
    score = cross_val_score(model, parsedData.x, y=parsedData.y,cv=10);
    print "   ",numpy.min(score)," < ",numpy.mean(score)," < ",numpy.max(score)

    # And the unaligned data
    parsedData = readData('UnalignedData.csv')
    model = svm.SVC(kernel='rbf',gamma=0.01, C=100)
    print "    10 Fold CVE C=100, Gamma = 0.01 ..."
    score = cross_val_score(model, parsedData.x, y=parsedData.y,cv=10);
    print "   ",numpy.min(score)," < ",numpy.mean(score)," < ",numpy.max(score)

    # And the boundary data
    parsedData = readData('BoundaryData.csv')
    model = svm.SVC(kernel='rbf',gamma=0.01, C=100)
    print "    10 Fold CVE C=100, Gamma = 0.01 ..."
    score = cross_val_score(model, parsedData.x, parsedData.y,cv=10);
    print "   ",numpy.min(score)," < ",numpy.mean(score)," < ",numpy.max(score)


    # And the advanced data
    parsedData = readData('BeamPatternData.csv')
    model = svm.SVC(kernel='rbf',gamma=0.01, C=100)
    print "    10 Fold CVE C=100, Gamma = 0.01 ..."
    score = cross_val_score(model, parsedData.x, parsedData.y,cv=10);
    print "   ",numpy.min(score)," < ",numpy.mean(score)," < ",numpy.max(score)

    # And the advanced data
    parsedData = readData('BeamShallowData.csv')
    model = svm.SVC(kernel='rbf',gamma=0.01, C=100)
    print "    10 Fold CVE C=100, Gamma = 0.01 ..."
    score = cross_val_score(model,parsedData.x,parsedData.y,cv=10);
    print "   ",numpy.min(score)," < ",numpy.mean(score)," < ",numpy.max(score)
    
