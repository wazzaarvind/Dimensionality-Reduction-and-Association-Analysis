

from __future__ import division
from collections import OrderedDict
from sklearn.manifold import TSNE
from sklearn.decomposition import TruncatedSVD
from sklearn.random_projection import sparse_random_matrix
import io
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def part1pca(fileName):	
	countTab=[]
	file=open(fileName,'r')
	line=file.readline()

	columnsLimit=len(line.strip().split("\t"))

	matrix=np.loadtxt(fileName,usecols=range(0,columnsLimit-1))
	
	mean=np.mean(matrix,axis=0)
	
	adjustedMatrix=np.subtract(matrix,mean)

	#Execute covariance manually as results of np.cov differ based pn python version
	covarianceMatrix=(1/matrix.shape[0])*np.dot(np.transpose(adjustedMatrix),adjustedMatrix) 
	
	eigenValues, eigenVectors=np.linalg.eig(covarianceMatrix)
	
	eigenValueFS=np.argsort(eigenValues,-1,kind='quicksort',order=None)[::-1][:2]
	
	eigenVectorsFS=eigenVectors[:,eigenValueFS]
	
	i=0
	arrayMultiple=np.empty([matrix.shape[0],2])
	for EV in eigenVectorsFS.T:
		arrayMultiple[:,i]=np.dot(adjustedMatrix,EV.T) 
		i=i+1

	diseasesArray=[]
	label_array=[]
	file2=open(fileName,'r')
	for line in file2:
		diseasesArray.append(line.strip().split("\t")[-1])

	#Define a dictionary which will hold unique integer as key for each unique disease as value
	dictionary = dict()
	counter=0 #This varaible will use 0 to n integers for each unique disease
	for disease in diseasesArray:
		if disease in dictionary:
			continue #Ignore disease if its already a part of the disease. We need only unique disease to assign unique integers
		dictionary[disease]=counter #If unique disease or first disease occurrance is found assign counter variable's current value against disease
		counter+=1 #Increament counter so as to pair next unique disease against an integer
		label_array.append(disease)

	#List of integers that holds the n different rows of diseases represented as a unique integer assigned for each unique disease
	numbers=[dictionary[disease] for disease in diseasesArray] 
	
	#Two eigen vectors are assigned to different variables and used to plot the scatter plot 
	x=arrayMultiple[:,0]
	y=arrayMultiple[:,1]
	
	#PCA Graph
	colors = cm.Set1(np.linspace(0, 1, len(label_array)))
	
	plt.figure()
	ax=plt.subplot(1,1,1)

	for l in range(0,len(label_array)):    # looping through class labels
		p1=[]	
		p2=[]
		for nums in range(0,len(x)):	# looping through points x 
			if numbers[nums]==l:        # if number belongs to the particular class label
			 	p1.append(x[nums])
			 	p2.append(y[nums])
		ax.scatter(p1,p2,s=50,color=colors[l] , label=label_array[l] )
					 	
	plt.title("Algorithm: PCA - "+"Dataset: "+fileName)	
	ax.legend(scatterpoints=1 , loc='best')
	plt.show()

	#Plotting graph for SVD FEATURE MATRIX
	U = TruncatedSVD(n_components=2).fit_transform(matrix)
	x2 = U[:,0]
	y2 = U[:,1]

	plt.figure()
	ax2=plt.subplot(1,1,1)
	for l in range(0,len(label_array)):    # looping through class labels
		p1=[]	
		p2=[]
		for nums in range(0,len(x2)):	# looping through points x 
			if numbers[nums]==l:        # if number belongs to the particular class label
			 	p1.append(x2[nums])
			 	p2.append(y2[nums])
		ax2.scatter(p1,p2,s=50,color=colors[l] , label=label_array[l] )
					 	
	plt.title("Feature Matrix; Algorithm: SVD - "+"Dataset: "+fileName)		
	ax2.legend(scatterpoints=1 , loc='best')
	plt.show()

	#Plotting graph for SVD MEAN CENTERED MATRIX 
	U = TruncatedSVD(n_components=2).fit_transform(adjustedMatrix)
	x3 = U[:,0]
	y3 = U[:,1]

	plt.figure()
	ax3=plt.subplot(1,1,1)
	for l in range(0,len(label_array)):    # looping through class labels
		p1=[]	
		p2=[]
		for nums in range(0,len(x3)):	# looping through points x 
			if numbers[nums]==l:        # if number belongs to the particular class label
			 	p1.append(x3[nums])
			 	p2.append(y3[nums])
		ax3.scatter(p1,p2,s=50,color=colors[l] , label=label_array[l] )
					 	
	plt.title("Mean Centered Matrix; Algorithm: SVD - "+"Dataset: "+fileName)		
	ax3.legend(scatterpoints=1 , loc='best')
	plt.show()

	#Plotting graph for t-SNE FEATURE MATRIX
	tsneResult = TSNE(n_components=2).fit_transform(matrix)
	x1 = tsneResult[:,0]
	y1 = tsneResult[:,1]

	plt.figure()
	ax1=plt.subplot(1,1,1)
	for l in range(0,len(label_array)):    # looping through class labels
		p1=[]	
		p2=[]
		for nums in range(0,len(x1)):	# looping through points x 
			if numbers[nums]==l:        # if number belongs to the particular class label
			 	p1.append(x1[nums])
			 	p2.append(y1[nums])
		ax1.scatter(p1,p2,s=50,color=colors[l] , label=label_array[l] )
					 	
	plt.title("Feature Matrix; Algorithm: TSNE - "+"Dataset: "+fileName)	
	ax1.legend(scatterpoints=1 , loc='best')
	plt.show()

	#Plotting graph for t-SNE MEAN CENTERED MATRIX 
	tsneResult = TSNE(n_components=2).fit_transform(adjustedMatrix)
	x4 = tsneResult[:,0]
	y4 = tsneResult[:,1]

	plt.figure()
	ax4=plt.subplot(1,1,1)
	for l in range(0,len(label_array)):    # looping through class labels
		p1=[]	
		p2=[]
		for nums in range(0,len(x4)):	# looping through points x 
			if numbers[nums]==l:        # if number belongs to the particular class label
			 	p1.append(x4[nums])
			 	p2.append(y4[nums])
		ax4.scatter(p1,p2,s=50,color=colors[l] , label=label_array[l] )
					 	
	plt.title("Mean Centered Matrix; Algorithm: TSNE - "+"Dataset: "+fileName)	
	ax4.legend(scatterpoints=1 , loc='best')
	plt.show()


#part1pca("pca_demo.txt")
part1pca("pca_a.txt")
part1pca("pca_b.txt")
part1pca("pca_c.txt")







