# Pattern-recognition
Implementation of different pattern recognition algorithms in MATLAB.

**Task 1 - Letter Classification**

The task requirement is to classify letters written on a graphics board, using a classifier 
based on hypothesis testing, and a parametric classifier. Firstly, we choose the features that we will use to recognize different classes of letters. For this, a special function, 
feature_extraction, is written, which finds 15 features from the information we have about each 
letter. The features we have singled out are min and max positions of shapes on x and y axes, 
ranges of these positions, min and max pressure, pressure range, mean speed value on the x 
and y axes, as well as the number of letters selected on the left, right upper and lower half of 
the image.
As it is very ungrateful to work with such a number of features, both because of the 
mathematical complexity, and because it is much easier for us to work with forms with 2 or 3 
features, because in this way we can visually display the shapes and more intuitively design 
the classifier, the next step is to use one of the methods of dimension reduction.
The first classifier we have designed is the Bayesian classifier, about which there will be more words in Task 2.
The last requirement of the task is to choose two letters that are linearly separable and choose one of the parametric classification methods for their classification. The letters chosen are E and S, and the parametric classification method is the wanted output method.
Please refer to Project report.pdf for more detailed explanation of the implementation of the solution of the task.

**Task 2 – Hypothesis Testing**

A special type of classifiers are classifiers formed based on hypothesis testing. In them, the 
classification is made by setting a hypothesis, which is then tested in the test, and based on the results of the set tests, we decide in which class to classify certain forms. Classifiers we desgined are the Bayesian classifier, the minimum price classifier, the Neyman-Pearson classifier and Wald's sequential classifier.
Advantages, disadvantages and implementation of all the classifiers are explaind in Project report.pdf.

**Task 3 – Parametric classification methods**

Since often the probability density functions of the classes are not available for classification of forms, the solution is found in new classifier design methods. One of these solutions are parametric classification methods, which do not require knowledge of the functions of the probability densities.
The first classifier we designed is a part-by-part linear classifier which we use in situations where we have inter-linearly separable classes, but their number is greater than 2, and one classification line is not enough. We design it by designing, in the first case, an optimal linear classifier between every two classes.
In the next step, the discriminatory lines of the part-by-part linear classifier are designed using the method of desired output.
The next task requirement is to generate two nonlinear separable classes. In order to successfully classify nonlinear separable classes, we need to switch from using a linear classifier to using a quadratic classifier. We designed the quadratic classifier by using the desired output method again, between every two classes.
Please refer to Project report.pdf for more detailed explanation of the implementation of the solution of the task.

**Task 4 - Clustering**

Clustering is a classification method in which we do not rely on any prior knowledge related 
to our samples. Unlike the previous method of classification, where we owned a training set on 
which we find the parameters of our model, we do not have a training set in clustering. 
Clustering methods can also be parametric and non-parametric, but only parametric methods 
were used in this task. The parametric methods used in the task are the c-mean method, the 
maximum credibility method, and the quadratic decomposition method.
Advantages, disadvantages and implementation of all the classifiers are explaind in Project report.pdf.
