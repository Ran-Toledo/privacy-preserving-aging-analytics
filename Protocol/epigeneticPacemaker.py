# import the data retrieval class
from EpigeneticPacemaker.ExampleData.DataSets import get_example_data
from typing import Dict, Tuple
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
from scipy import stats
import math
# from scipy.stats import chi2_contingency
import subprocess


# clear from data variables correspond to small ages
def clear_data_by_ages(ages, methylation_values, limit):
    x = 0
    while x < len(ages):
        if ages[x] < limit:
            # delete age and site values related to it
            methylation_values = np.delete(methylation_values, x, axis=1)
            ages = np.delete(ages, x)
        else:
            x += 1
    return ages, methylation_values


# func for calculating absolute Pearson correlation
def pearson_correlation(ages, methylation_values):
    p_values = []
    for x in range(len(methylation_values)):
        p_values.append(abs(scipy.stats.pearsonr(methylation_values[x], ages)[0]))
    return p_values


# create y vector for sites with high Pearson correaltion
def create_y_vector(s):
    y = []
    for i in range(len(s)):
        for j in range(len(s[i])):
            y.append(s[i][j])  # entry Si,j
    return y


# Naive func for calculating rates and start states(βˆ = (X^T *X)^−1 * X^T*y)
# gets : ages (X matrix), sites (y vecrtor), n (the number of sites)
def calculate_rates_and_start_states_naive(ages, sites, n):
    m = len(ages)
    matrix = []
    # create matrix X
    for x in range(n):
        for y in range(m):
            row = []
            place = 0
            while (place < x):
                row.append(0)
                place += 1
            row.append(ages[y])
            place += 1
            while place < n:
                row.append(0)
                place += 1
            place = 0
            while (place < x):
                row.append(0)
                place += 1
            row.append(1)
            place += 1
            while place < n:
                row.append(0)
                place += 1
            matrix.append(row)
    mat = []
    for x in range(len(matrix)):
        temp = []
        for y in range(len(matrix[0])):
            temp.append(matrix[x][y])
        mat.append(temp)
    f = open("demofile2.txt", "w")
    f.write(str(mat))
    f.close()
    # create X^T
    matrix_transpose = np.array(matrix).transpose()
    # multipe X^T *X
    mul = np.matmul(matrix_transpose, np.array(matrix))
    # get matrix inverse ((X^T *X)^−1)
    inverse = np.linalg.inv(mul)
    # multiple [(X^T *X)^−1] and  [X^T]
    mul = np.matmul(inverse, matrix_transpose)
    # multiple the result by vector y
    result = np.matmul(mul, np.array(sites))
    ri = []
    s0i = []
    # get values for ri and si0
    for x in range(int(len(result) / 2)):
        ri.append(result[x])
    for x in range(int(len(result) / 2), len(result)):
        s0i.append(result[x])
    return ri, s0i


# advanced func for calculating rates and start states(βˆ = (X^T *X)^−1 * X^T*y)
# gets : ages and sites
def calculate_rates_and_start_states(ages, sites):
    a = sum(ages)
    b = 0
    for i in range(len(ages)):
        b += ages[i] * ages[i]
    temp = a * a
    Lambda = 1 / (temp - len(ages) * b)
    first_vector = []
    second_vector = []
    for i in range(len(ages)):
        temp = Lambda * (-len(ages) * ages[i] + a)
        first_vector.append(temp)
    for i in range(len(ages)):
        temp = Lambda * (a * ages[i] - b)
        second_vector.append(temp)
    ri = []
    i = 0
    while i < len(sites):
        multiple_vectors = 0
        for j in range(len(first_vector)):
            multiple_vectors += sites[i] * first_vector[j]
            i += 1
        ri.append(multiple_vectors)
    i = 0
    s0i = []
    while i < len(sites):
        multiple_vectors = 0
        for j in range(len(second_vector)):
            multiple_vectors += sites[i] * second_vector[j]
            i += 1
        s0i.append(multiple_vectors)
    return ri, s0i


# calculate epigenetic ages by tj = ∑(i≤n)[ri*(si,j-s0i))]/∑(i≤n)[ri^2]
def calculate_ages(S, s0i, ri):
    denominator = 0
    for i in range(len(ri)):
        denominator += ri[i] * ri[i]
    tj = []
    for j in range(len(S[0])):
        val = 0
        for i in range(len(S)):
            val += ri[i] * (S[i][j] - s0i[i])
        tj.append(val / denominator)
    return tj


# calculate RSS,  (RSS =∑(i≤n)∑(j≤m)((Si,j – (s0i + ritj))^2).
def calculate_rss(s, s0i, ri, tj):
    rss = 0
    for i in range(len(ri)):
        for j in range(len(tj)):
            temp = s[i][j] - (s0i[i] + ri[i] * tj[j])
            rss += temp * temp
    return rss


# EPM algorithm
# gets S' matrix, Y vector, chronological ages - tj, delta = minimal improvement acceptable,  n = maximum number of iterations
def EPM(S, y, tj, delta, n):
    m = 2
    # first iteration
    ri, s0i = calculate_rates_and_start_states(tj, y)
    tj = calculate_ages(S, s0i, ri)
    RSS0 = calculate_rss(S, s0i, ri, tj)
    if (n == 1):
        return ri, s0i, tj, RSS0
    # second iteration
    ri, s0i = calculate_rates_and_start_states(tj, y)
    RSS1 = calculate_rss(S, s0i, ri, tj)
    tj = calculate_ages(S, s0i, ri)
    # rest of iterations
    while (RSS0 - RSS1 > delta and m < n):
        RSS0 = RSS1
        ri, s0i = calculate_rates_and_start_states(tj, y)
        tj = calculate_ages(S, s0i, ri)
        RSS1 = calculate_rss(S, s0i, ri, tj)
        m += 1
    return ri, s0i, tj, RSS1, m


# create scatter graph
def create_graph(x, y):
    plt.scatter(x, y, color="blue", s=10)
    xval = min(x)
    yval = min(y)
    # xticks with jumps of 5, add 0 val, max and min vals
    xaxes = []
    # yticks with jumps of 5, add 0 val, max and min vals
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
    plt.show()


# calculate chi square and p value
def chi_square(num_individuals, num_sites, val1, val2):
    # calc chi square
    chi = num_individuals * num_sites * math.log(val1 / val2)
    # calc p value, use chi value and degree of freedom
    p_val = 1 - stats.chi2.cdf(chi, num_individuals)
    return chi, p_val


# calculate noise, equ : sij - (s0i + ritj)  = εi,j
def noise(S, s0, r, t):
    total = 0
    for i in range(len(S)):
        for j in range(len(S[i])):
            noise = abs(S[i][j] - s0[i] - r[i] * t[j])
            total += noise
    return total / (len(S) * len(S[0]))


# ratio between its MC to PM rate
def graph_ratio_sites(rMC, rPM):
    arr = []
    s = 0
    for x in range(len(rMC)):
        arr.append(rMC[x] / rPM[x])
        if rMC[x] < rPM[x]:
            s += 1
    arr.sort()
    plt.plot(arr)
    plt.axhline(1, color='red')
    plt.show()


# ratio between epigenetic to chronological age
def graph_ratio_ages(e_age, c_age):
    arr = []
    x = []
    for y in range(len(e_age)):
        arr.append(c_age[y] / e_age[y])
        x.append(y)
    sort = [x for _, x in sorted(zip(c_age, arr))]
    plt.plot(x, sort)
    plt.gca().yaxis.grid(True)
    plt.show()


if __name__ == '__main__':
    # retrieve the training and testing data
    data1, data2 = get_example_data()
    # unpack the training and testing data
    samples1, cpg_sites1, ages1, methylation_values1 = data1
    samples2, cpg_sites2, ages2, methylation_values2 = data2
    ages = np.hstack([ages1, ages2])
    methylation_values = np.hstack([methylation_values1, methylation_values2])
    # take only peaple with age >= 0
    # ages, methylation_values = Clear_data_by_ages (ages, methylation_values,0)
    # get the absolute value of the correlation coefficient
    pearson = np.array(pearson_correlation(ages, methylation_values))
    # return list of site indices with a high absolute correlation coefficient
    sites = np.where(pearson > .85)[0]
    print("the number of individuals is : ", str(len(ages)))
    print("the number of sites is : ", str(len(sites)))
    # creat S matrix :
    S = methylation_values[sites]
    # creat y vector :
    y = create_y_vector(S)
    # run EPM algorithm
    ri, s0i, tj, RSS, m = EPM(S, y, ages, 0.00001, 1000)
    print("number of iterrations is: ", str(m))
    create_graph(ages, tj)


