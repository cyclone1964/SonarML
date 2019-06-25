# This is an attempt to apply an SVM to the data. Specficially, we
# train 4 SVMS each 48 bins wide plus annotations.
import numpy as numpy
from sklearn import svm
from sklearn.utils import resample
from sklearn.metrics import accuracy_score
from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import RandomForestClassifier

# Define the file reader
execfile('readData.py')
    
# Read in the data, buld the annotations, and do the normalization
allData = readData('BoundaryData.csv')
y = allData.y
x = numpy.concatenate((allData.annotations,allData.x),axis=1);
for index in range(len(y)):
    limits = [numpy.min(x[index,4:]),
              numpy.max(x[index,4:])]
    x[index,4:] = (x[index,4:] - limits[0])/(limits[1]-limits[0])

numTraining = 6000
svmScores = []
forestScores = []
annScores = []

# Set up a list of indices
indices = [i for i in range(len(y))]
    
svm = svm.SVC(kernel='rbf',gamma=0.001, C=100)
forest = RandomForestClassifier()
ann = MLPClassifier(activation='logistic',
                      solver='lbfgs',
                      hidden_layer_sizes = (64))


for iterationCount in range(200):
    print "Iteration: ",iterationCount,

    # Get the training indices by resampling with replacement, and
    # then the test indices are the ones that are left. I got this
    # code off the net I admit it
    trainingIndices = resample(indices,n_samples = numTraining)
    testIndices = numpy.array([i for i in indices if i not in trainingIndices])

    svm.fit(x[trainingIndices,:],y=y[trainingIndices])
    predictions = svm.predict(x[testIndices,:])
    score = accuracy_score(y[testIndices],predictions)
    svmScores.append(score)
    print "svmScore ",score,

    forest.fit(x[trainingIndices,:],y=y[trainingIndices])
    predictions = forest.predict(x[testIndices,:])
    score = accuracy_score(y[testIndices],predictions)
    forestScores.append(score)
    print "forestScore ",score,

    ann.fit(x[trainingIndices,:],y=y[trainingIndices])
    predictions = ann.predict(x[testIndices,:])
    score = accuracy_score(y[testIndices],predictions)
    annScores.append(score)
    print "annScore ",score
    


svmScores.sort()
forestScores.sort()
annScores.sort()

print "svm Confidence [",svmScores[10]," - ", svmScores[190]
print "forest Confidence [",forestScores[10]," - ", forestScores[190]
print "ann Confidence [",annScores[10]," - ", annScores[190]
    
