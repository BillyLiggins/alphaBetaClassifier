"""
You have found that this simple logistic regression performs a lot better when you train on the posr and BAB classifier as oppose to the energy. 

This is werid!

however may be explained by the .... Think about it. 

"""


import numpy as np
import matplotlib.pyplot as plt
from sklearn import linear_model, datasets

data_bi=np.loadtxt("Classifier_data_bi.dat",skiprows=1,delimiter=",",usecols=(0,2))

target_bi=np.ones(data_bi.shape[0])

data_po=np.loadtxt("Classifier_data_po.dat",skiprows=1,delimiter=",",usecols=(0,2))

target_po=np.zeros(data_po.shape[0])

data=np.array([y for x in map(None,data_bi,data_po[:data_bi.shape[0]]) for y in x if y is not None])

target=np.array([y for x in map(None,target_bi,target_po[:data_bi.shape[0]]) for y in x if y is not None]) 

print "Data loaded"

X= data[:data.shape[0]/2]
Y= target[:target.shape[0]/2]

testData= data[data.shape[0]/2:]
testTarget= target[target.shape[0]/2:]


h = 2  # step size in the mesh

logreg = linear_model.LogisticRegression(C=1e5)

# we create an instance of Neighbours Classifier and fit the data.
logreg.fit(X, Y)

# Plot the decision boundary. For that, we will assign a color to each
# point in the mesh [x_min, m_max]x[y_min, y_max].
x_min, x_max = X[:, 0].min() - .5, X[:, 0].max() + .5
y_min, y_max = X[:, 1].min() - .5, X[:, 1].max() + .5
xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))
ZZ = logreg.predict(np.c_[xx.ravel(), yy.ravel()])
Z = logreg.predict(testData)

print Z

data_plot=(testTarget==Z).astype(int)
print data_plot

# plt.plot(data_plot)
# plt.savefig("compareResults.png")
success=data_plot.tolist().count(1)
fail=data_plot.tolist().count(0)
print "Number of successes = "+ str(success)
print "Number of failures = "+ str(fail)
print "chance of failure (fail/total) = " +str(float(fail)/(fail+success))



# Put the result into a color plot
ZZ = ZZ.reshape(xx.shape)
plt.figure(1, figsize=(4, 3))
plt.pcolormesh(xx, yy, ZZ, cmap=plt.cm.Paired)

#Plot also the training points
plt.scatter(X[:, 0], X[:, 1], c=Y, edgecolors='k', cmap=plt.cm.Paired)
plt.xlabel('MC Energy (MeV)')
plt.ylabel('BAB')

plt.xlim(xx.min(), xx.max())
plt.ylim(yy.min(), yy.max())
plt.xticks(())
plt.yticks(())

plt.show()
