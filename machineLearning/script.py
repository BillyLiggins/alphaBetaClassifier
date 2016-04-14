import numpy as np
import matplotlib.pyplot as plt
from sklearn import linear_model, datasets

data_bi=np.loadtxt("Classifier_data_bi.dat",skiprows=1,delimiter=",")

target_bi=np.ones(data_bi.shape[0])

data_po=np.loadtxt("Classifier_data_po.dat",skiprows=1,delimiter=",")

target_po=np.zeros(data_po.shape[0])

data=np.array([y for x in map(None,data_bi,data_po[:data_bi.shape[0]]) for y in x if y is not None])#

target=np.array([y for x in map(None,target_bi,target_po[:data_bi.shape[0]]) for y in x if y is not None]) 

X=trainingData= data[:data.shape[0]/2]
Y=trainingTarget= target[:target.shape[0]/2]

testData= data[data.shape[0]/2:]
testTarget= target[target.shape[0]/2:]

# import some data to play with
# iris = datasets.load_iris()
# X = iris.data[:, :2]  # we only take the first two features.
# Y = iris.target

h = .02  # step size in the mesh

logreg = linear_model.LogisticRegression(C=1e5)

# we create an instance of Neighbours Classifier and fit the data.
logreg.fit(X, Y)

# Plot the decision boundary. For that, we will assign a color to each
# point in the mesh [x_min, m_max]x[y_min, y_max].
x_min, x_max = X[:, 0].min() - .5, X[:, 0].max() + .5
y_min, y_max = X[:, 1].min() - .5, X[:, 1].max() + .5
xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))
Z = logreg.predict(np.c_[xx.ravel(), yy.ravel()])

# Put the result into a color plot
Z = Z.reshape(xx.shape)
plt.figure(1, figsize=(4, 3))
plt.pcolormesh(xx, yy, Z, cmap=plt.cm.Paired)

# Plot also the training points
plt.scatter(X[:, 1], X[:, 2], c=Y, edgecolors='k', cmap=plt.cm.Paired)
plt.xlabel('MC Energy (MeV)')
plt.ylabel('BAB')

plt.xlim(xx.min(), xx.max())
plt.ylim(yy.min(), yy.max())
plt.xticks(())
plt.yticks(())

plt.show()
