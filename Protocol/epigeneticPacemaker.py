# import the data retrieval class
from EpigeneticPacemaker.ExampleData.DataSets import get_example_data
from typing import Dict, Tuple
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
from scipy import stats
import math
#from scipy.stats import chi2_contingency
import subprocess
import numpy as np

#clear from data variables correspond to small ages 
def Clear_data_by_ages (test_ages, test_methylation_values, limit):
    x = 0
    while x < len(test_ages):
        if (test_ages[x] < limit):
            #delete age and site values related to it
            test_methylation_values = np.delete(test_methylation_values, x, axis=1)
            test_ages = np.delete(test_ages, x)
        else:
             x+=1
    return test_ages, test_methylation_values

#----------------------------------------------------------------------------------------
#func for calculating absolute Pearson correaltion 
def pearson_correlation (test_ages, test_methylation_values):
    p_values = []
    for x in range (len(test_methylation_values)):
       p_values.append(abs(scipy.stats.pearsonr(test_methylation_values[x], test_ages)[0])) 
    return p_values
    
#create y vector for sites with high Pearson correaltion
def Create_Y_vector (S):
     y = []
     for i in range (len(S)):
        for j in range (len(S[i])):
            y.append(S[i][j])     #entry Si,j
     return y
    
#create S' entries for sites with high Pearson correaltion    
def Create_S_matrix (training_sites, test_methylation_values):
     sites = []
     for x in range (len(training_sites)):
            sites.append(test_methylation_values[training_sites[x]])     #entry Si,j
     return sites
    
#Naive func for calculating rates and start states(βˆ = (X^T *X)^−1 * X^T*y)
#gets : ages (X matrix), sites (y vecrtor), n (the number of sites)
def calculate_rates_and_startStates_Naive (ages, sites,n) :
    m = len(ages)
    matrix = []
    #create matrix X
    for x in range (n):
        for y in range (m):
            row = []
            place = 0
            while (place < x):
                row.append(0)
                place+=1
            row.append(ages[y])
            place+=1
            while place<n:
              row.append(0)
              place+=1
            place = 0
            while (place < x):
                row.append(0)
                place+=1
            row.append(1)
            place+=1
            while place<n:
              row.append(0)
              place+=1        
            matrix.append(row)
    mat = []
    for x in range (len(matrix)):
        temp = []
        for y in range (len(matrix[0])):
                            temp.append(matrix[x][y])
        mat.append(temp)    
    #f = open("demofile2.txt", "w")
    #f.write(str(mat))
    #f.close()
    #create X^T       
    matrix_transpose = np.array(matrix).transpose()
    #multipe X^T *X
    mul =np.matmul(matrix_transpose, np.array(matrix))
    #get matrix inverse ((X^T *X)^−1)
    inverse = np.linalg.inv(mul)
    #multiple [(X^T *X)^−1] and  [X^T]
    mul =np.matmul(inverse,matrix_transpose)
    #multiple the result by vector y
    result = np.matmul(mul, np.array(sites))
    ri = []
    s0i = []
    #get values for ri and si0
    for x in range (int(len(result)/2)):
        ri.append(result[x])
    for x in range (int(len(result)/2), len(result)):    
        s0i.append(result[x])
    return ri, s0i

def calculate_rates_and_startStates(ages, sites) :
    a = sum(ages)
    b = 0
    for i in range (len(ages)):
        b += ages[i] * ages[i]    
    temp = a * a
    Lambda = 1 / (temp - len(ages)*b)
    first_vector = []
    second_vector = []
    for i in range (len(ages)):
       temp = Lambda * (-len(ages)*ages[i] + a)
       first_vector.append(temp)
    for i in range (len(ages)):
       temp =  Lambda* (a*ages[i] - b)
       second_vector.append(temp)
    ri = []
    i = 0
    while i < len(sites):
        multiple_vectors = 0
        for j in range (len(first_vector)):
            multiple_vectors += sites[i]*first_vector[j]
            i+=1
        ri.append( multiple_vectors)
    i = 0
    s0i = []
    while i < len(sites):
        multiple_vectors = 0
        for j in range (len(second_vector)):
            multiple_vectors += sites[i]*second_vector[j]
            i+=1
        s0i.append( multiple_vectors)
    return ri, s0i
# calculate epigenetic ages by tj = ∑(i≤n)[ri*(si,j-s0i))]/∑(i≤n)[ri^2]
def calcultae_eAges (S, s0i, ri):
    denominator = 0
    for i in range (len(ri)):
        denominator += ri[i] * ri[i]
    tj = []
    for j in range (len(S[0])):
        val = 0
        for i in range (len(S)):
            val += ri[i] * (S[i][j] - s0i[i])
        tj.append(val/denominator)
    return tj    
    
#calculate RSS,  (RSS =∑(i≤n)∑(j≤m)((Si,j – (s0i + ritj))^2). 
def calculateRSS (S, s0i, ri, tj):
    RSS = 0
    for i in range (len(ri)):
        for j in range (len(tj)):
            temp = S[i][j] - (s0i[i] + ri[i]*tj[j])
            RSS += temp*temp
    return RSS

#EPM algorithm
#gets S' matrix, Y vector, chronological ages - tj, delta = minimal improvement acceptable,  n = maximum number of iterations
def EPM (S, y, tj, delta, n):
    m = 2
    #first iteration 
    ri, s0i = calculate_rates_and_startStates_Naive(tj, y, len(S))
    #print (ri)
    #print (s0i)
    tj = calcultae_eAges(S, s0i, ri)
    RSS0 = calculateRSS (S, s0i, ri, tj)
    #second iteration 
    ri, s0i = calculate_rates_and_startStates_Naive(tj, y, len(S))
    tj = calcultae_eAges(S, s0i, ri)
    RSS1 = calculateRSS (S, s0i, ri, tj)
    #rest of iterations
    while (RSS0 - RSS1 > delta and m < n):

  #  while (m < n):
        RSS0 = RSS1
        ri, s0i = calculate_rates_and_startStates_Naive(tj, y, len(S))
        tj = calcultae_eAges(S, s0i, ri)
        RSS1 = calculateRSS (S, s0i, ri, tj)
        m+=1
    #print (m)
    return ri, s0i,tj
'''
#create graphs and calculate R-squared 
def create_graph(x,y):
  plt.scatter(x,y, color = "red", s = 15)
  xval = min(x)
  yval = min(y)
  #xticks with jumps of 5, add 0 val, max and min vals
  xaxes = []
  #yticks with jumps of 5, add 0 val, max and min vals
  yaxes = []
  while (xval < max(x)):
    xaxes.append(xval)
    xval += 5
  if (0 not in xaxes):
      xaxes.append(0)
  xaxes.append(max(x))
  while (yval < max(y)):
    yaxes.append(yval)
    yval += 5
  if (0 not in yaxes):
      yaxes.append(0)
  yaxes.append(max(y)) 
  plt.yticks(yaxes)
  plt.xticks(xaxes, rotation='vertical')
  plt.xlabel("chronological ages")
  plt.ylabel("epigenetic ages")
  plt.gca().xaxis.grid(True)
  plt.gca().yaxis.grid(True)
  res = stats.linregress(x, y)
  print(f"R-squared: {res.rvalue**2:.6f}") 
  plt.plot(np.array(x), res.intercept + res.slope*np.array(x), 'b', label=f"R-squared: {res.rvalue**2:.6f}")
  z = np.polyfit(x.flatten(), np.array(y).flatten(), 1)
  plt.title("y=%.6fx+%.6f"%(z[0],z[1])) 
  plt.legend()
  plt.show()
'''
def create_graph(x, y):
    #logging.disable(logging.DEBUG)
    plt.scatter(x, y, color="green", s=15)
    plt.xlabel("Chronological ages")
    plt.ylabel("Epigenetic ages")
    plt.gca().xaxis.grid(True)
    plt.gca().yaxis.grid(True)
    res = stats.linregress(x, y)
    plt.plot(np.array(x), res.intercept + res.slope * np.array(x), 'b', label=f"R-squared: {res.rvalue ** 2:.6f}")
    # z = np.polyfit(np.array(x).flatten(), np.array(y).flatten(), 1)
    # plt.title("y=%.6fx+%.6f"%(z[0],z[1]))
    plt.title("Epigenetic PaceMaker")
    plt.legend()
    plt.show()
def graph (x,y):
  plt.scatter(x,y, color = "red", s = 15)
  xval = min(x)
  yval = min(y)
  #xticks with jumps of 5, add 0 val, max and min vals
  xaxes = []
  #yticks with jumps of 5, add 0 val, max and min vals
  yaxes = []
  while (xval < max(x)):
    xaxes.append(xval)
    xval += 5
  if (0 not in xaxes):
      xaxes.append(0)
  xaxes.append(max(x))
  while (yval < max(y)):
    yaxes.append(yval)
    yval += 5
  if (0 not in yaxes):
      yaxes.append(0)
  yaxes.append(max(y)) 
  yaxes.append(max(x)) 
  plt.yticks(yaxes)
  plt.xticks(xaxes, rotation='vertical')
  plt.xlabel("chronological ages")
  plt.ylabel("epigenetic ages")
  plt.gca().xaxis.grid(True)
  plt.gca().yaxis.grid(True)
  x = np.array(x) 
  y = np.array(y)
  n = np.size(x)  
  x_mean = np.mean(x)
  y_mean = np.mean(y)
  Sxy = np.sum(x*y)- n*x_mean*y_mean
  Sxx = np.sum(x*x)-n*x_mean*x_mean  
  b1 = Sxy/Sxx
  b0 = y_mean-b1*x_mean
  print('slope b1 is', b1)
  print('intercept b0 is', b0)
  y_pred = b1 * x + b0
  plt.plot(x, y_pred, color = 'green')
  res = stats.linregress(x, y)
  print(f"R-squared: {res.rvalue**2:.6f}")
  plt.plot(np.array(x), res.intercept + res.slope*np.array(x), 'b', label=f"R-squared: {res.rvalue**2:.6f}")
  plt.show()
def main ():
    
    # retrieve the training and testing data
    test_data, train_data = get_example_data()
    # unpack the training and testing data
    test_samples, test_cpg_sites, test_ages, test_methylation_values = test_data
    train_samples, train_cpg_sites, train_ages, train_methylation_values = train_data
    # get the absolute value of the correlation coefficient
    pearson = np.array(pearson_correlation (train_ages, train_methylation_values))
    # return list of site indices with a high absolute correlation coefficient
    training_sites = np.where(pearson > .90)[0]
    training_sites = training_sites[:5]
    train_ages = train_ages[:24]
    #creat S matrix : 
    S = train_methylation_values[training_sites, :24]
    #creat S matrix : 
    # S = Create_S_matrix(training_sites, train_methylation_values)
    #creat y vector : 
    y = Create_Y_vector(S)
    #print (y)
    #run EPM algorithm
    #print(calculate_rates_and_startStates_Naive(train_ages, y, len(S)))
    ri, s0i,tj= EPM(S,y,train_ages, 0.00001, 1000)
    #print (tj)
    create_graph(train_ages, tj)
    return tj
    #run MC
    #riMC, s0iMC = calculate_rates_and_startStates(test_ages, y)
    #riMCc, s0iMCc =calculate_rates_and_startStates(test_ages, y)
    #RSSMC = calculateRSS(S, s0iMC, riMC, test_ages)
    #RSSEPM = calculateRSS(S, s0i, ri, tj)
    #calcultate chi squared
    #chi = len(test_ages) * len(test_methylation_values) * math.log(RSSMC/RSSEPM)
    #print(chi)
    #creat S matrix for train data : 
    #S = Create_S_matrix(training_sites, train_methylation_values)
    #calculate e-ages for train data
    #tj = calcultae_eAges(S, s0i, ri)
    #create graph and calculate r 
    #create_graph(train_ages, tj)

main()
