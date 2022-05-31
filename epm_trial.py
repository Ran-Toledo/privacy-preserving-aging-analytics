# Load libraries
from pandas import read_csv
from pandas.plotting import scatter_matrix
from matplotlib import pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC

import math as md
import statistics as sd

# import the data retrieval class
from EpigeneticPacemaker.ExampleData.DataSets import get_example_data
# retrieve the training and testing data
test_data, train_data = get_example_data()
# unpack the training and testing data
test_samples, test_cpg_sites, test_ages, test_methylation_values = test_data
train_samples, train_cpg_sites, train_ages, train_methylation_values = train_data




import numpy as np

def pearson_correlation(meth_matrix: np.array, phenotype: np.array) -> np.array:
    """calculate pearson correlation coefficient between rows of input matrix and phenotype"""
    # calculate mean for each row and phenotype mean
    matrix_means = np.mean(meth_matrix, axis=1)
    phenotype_mean = np.mean(phenotype)

    # subtract means from observed values
    transformed_matrix = meth_matrix - matrix_means.reshape([-1,1])
    transformed_phenotype = phenotype - phenotype_mean

    # calculate covariance
    covariance = np.sum(transformed_matrix * transformed_phenotype, axis=1)
    variance_meth = np.sqrt(np.sum(transformed_matrix ** 2, axis=1))
    variance_phenotype = np.sqrt(np.sum(transformed_phenotype ** 2))

    return covariance / (variance_meth * variance_phenotype)

# get the absolute value of the correlation coefficient
abs_pcc_coefficients = abs(pearson_correlation(train_methylation_values, train_ages)) 

# return list of site indices with a high absolute correlation coefficient
training_sites = np.where(abs_pcc_coefficients > .85)[0]





import matplotlib.pyplot as plt
from matplotlib import rc
from scipy import optimize
import scipy.stats as stats

# use latex formatting for plots
rc('text', usetex=True)

def r2(x,y):
    # return r squared
    return stats.pearsonr(x,y)[0] **2

def plot_known_predicted_ages(known_ages, predicted_ages, label=None):
    # define optimization function
    def func(x, a, b, c):
        return a * np.asarray(x)**0.5 + c
    # fit trend line
    popt, pcov = optimize.curve_fit(func, [1 + x for x in known_ages], predicted_ages)
    # get r squared
    rsquared = r2(predicted_ages, func([1 + x for x in known_ages], *popt))
    # format plot label
    plot_label = f'$f(x)={popt[0]:.2f}x^{{1/2}} {popt[2]:.2f}, R^{{2}}={rsquared:.2f}$'
    # initialize plt plot
    fig, ax = plt.subplots(figsize=(12,12))
    # plot trend line
    ax.plot(sorted(known_ages), func(sorted([1 + x for x in known_ages]), *popt), 'r--', label=plot_label)
    # scatter plot
    ax.scatter(known_ages, predicted_ages, marker='o', alpha=0.8, color='k')
    ax.set_title(label, fontsize=18)
    ax.set_xlabel('Chronological Age', fontsize=16)
    ax.set_ylabel('Epigenetic State', fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.legend(fontsize=16)
    plt.show()



from EpigeneticPacemaker.EpigeneticPacemaker import EpigeneticPacemaker

# initialize the EPM model 
epm = EpigeneticPacemaker(iter_limit=100, error_tolerance=0.00001)

# fit the model using the training data
epm.fit(train_methylation_values[training_sites,:], train_ages)

# generate predicted ages using the test data
test_predict = epm.predict(test_methylation_values[training_sites,:])

# plot the model results 
plot_known_predicted_ages(test_ages, test_predict, 'EpigeneticPacemaker')

