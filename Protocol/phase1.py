# # -*- coding: utf-8 -*-
# """
# Created on Sat May 28 13:20:04 2022
#
# @author: rtoledo
# """
#
# from phe import paillier
# import pandas as pd
# import numpy as np
# import math as md
#
#
# from EpigeneticPacemaker.ExampleData.DataSets import get_example_data
# from EpigeneticPacemaker.EpigeneticPacemaker import EpigeneticPacemaker
#
#
# def pearson_correlation(meth_matrix: np.array, phenotype: np.array) -> np.array:
#     """calculate pearson correlation coefficient between rows of input matrix and phenotype"""
#     # calculate mean for each row and phenotype mean
#     matrix_means = np.mean(meth_matrix, axis=1)
#     phenotype_mean = np.mean(phenotype)
#
#     # subtract means from observed values
#     transformed_matrix = meth_matrix - matrix_means.reshape([-1,1])
#     transformed_phenotype = phenotype - phenotype_mean
#
#     # calculate covariance
#     covariance = np.sum(transformed_matrix * transformed_phenotype, axis=1)
#     variance_meth = np.sqrt(np.sum(transformed_matrix ** 2, axis=1))
#     variance_phenotype = np.sqrt(np.sum(transformed_phenotype ** 2))
#
#     return covariance / (variance_meth * variance_phenotype)
#
#
# def Split_Data_to_different_owners(m):
#     # retrieve the training and testing data
#     test_data, train_data = get_example_data()
#     # unpack the training and testing data
#     #test_samples, test_cpg_sites, test_ages, test_methylation_values = test_data
#     train_samples, train_cpg_sites, train_ages, train_methylation_values = train_data
#
#     # get the absolute value of the correlation coefficient
#     abs_pcc_coefficients = abs(pearson_correlation(train_methylation_values, train_ages))
#
#     # return list of site indices with a high absolute correlation coefficient
#     training_sites = np.where(abs_pcc_coefficients > .85)[0]
#
#     train_methylation_values = train_methylation_values[training_sites,:]
#     '''
#     #Rounding up the data
#     print(len(train_methylation_values))
#
#     for i in range(0,len(train_ages)):
#         train_ages[i] = round(train_ages[i],2)
#
#     for i in range(0, len(train_methylation_values)):
#         for j in range(0, len(train_methylation_values[i])):
#             train_methylation_values[i][j] = round(train_methylation_values[i][j],2)
#
#     print("***Rounding up data complete.")
#     '''
#     met_values_split = np.hsplit(train_methylation_values, m)
#     train_ages = np.array(train_ages)
#     ages_split = np.split(train_ages, m)
#
#     train_samples = np.array(train_samples)
#     samples_split = np.split(train_samples, m)
#
#     data_owners = []
#     for i in range(0, m):
#         data_owners.append([samples_split[i].tolist(), train_cpg_sites, ages_split[i], met_values_split[i]])
#
#     print("***Splitting data complete.")
#
#     return data_owners
#
#
# def Encrypt_Data(DO_i, public_key):
#     train_samples, train_cpg_sites, train_ages, train_methylation_values = DO_i
#     cipher_ages = []
#     cipher_meth_values = []
#
#     for gsm_age in train_ages:
#         cipher_ages.append(public_key.encrypt(gsm_age))
#
#     print("Encrypting ages complete.")
#
#     for site in train_methylation_values:
#         cipher_site = []
#         for gsm_site_value in site:
#             cipher_site.append(public_key.encrypt(gsm_site_value))
#         cipher_meth_values.append(cipher_site)
#
#     print("Encrypting methylation site rates complete.")
#
#     return [cipher_ages, cipher_meth_values]
#
#
# class CSP:
#     def __init__(self):
#         self.public_key, self.__secret_key = paillier.generate_paillier_keypair()
#
#     def get_pk(self):
#         return self.public_key
#
#     '''
#     def decrypt(self, num):
#         return self.__secret_key.decrypt(num)
#     '''
#
#
# class MLE:
#     encrypted_data=[]
#     X=[]
#     beta=[]
#     y=[]
#     def __init__(self,pk):
#         self.pk=pk
#
#     def Send_to_MLE(self, encrypted_data):
#         self.encrypted_data.append(encrypted_data)
#
#     def getXandY(self):
#         """
#         Reminder:   X = ages matrix
#                     y = site methylation rates
#                     Model: (X^T*X)^-1 * X^T * y = beta
#         """
#         #Unifying the data
#         all_ages = []
#         all_meth_values = []
#         for ages, meth_values in encrypted_data:
#             all_ages.append(ages)
#             all_meth_values.append(meth_values)
#
#         #Defining X matrix
#         all_meth_values = np.array([np.array(meth_values) for meth_values in all_meth_values])
#         m, n = all_meth_values.shape
#         for i in range(0, m*n):
#             age_index = 0
#             X_row = []
#             for j in range(0, 2*n):
#                 if i == j:
#                     X_row.append(all_ages[age_index].ciphertext)
#                     age_index += 1
#                     if age_index==len(all_ages):
#                         age_index = 0
#                 if i+n == j:
#                     X_row.append(1)
#                 else:
#                     X_row.append(0)
#             self.X.append(X_row)
#
#         #Defining y vector
#         #y = [S_1,1 ... S_n,m] (iterating sites first)
#         transpose = all_meth_values.transpose()
#         for meth_site in transpose:
#             for gsm_value in meth_site:
#                 self.y.append(gsm_value)
#
#
#
# """
# Phase1:
#     1.      CSP generates the key pair (sk,pk), stores sk and makes pk public
# 	2.      Each data owner DO_i sends MLE ciphertexts computed using pk and the values in D_i
# 	3.      MLE uses these ciphertexts and the homomorphic property of the encryption scheme
#             in order to obtain encryptions of A and b (coefficient matrix and vector)
#             A = X   ,   b = y
# """
# #Step 1
# MyCSP = CSP()
# public_key = MyCSP.get_pk()
#
# #Step 2
# MyMLE = MLE(public_key)
# data_owners = Split_Data_to_different_owners(8)
#
# for DO_i in data_owners:
#     print("Encrypting Data Owner #" + (str)(data_owners.index(DO_i)+1) + "'s  data")
#     encrypted_data = Encrypt_Data(DO_i, public_key)
#     MyMLE.Send_to_MLE(encrypted_data)
#
# #Step 3
# MyMLE.getXandY()
#
#
#
#
# """
# Phase2:
#     1.     	MLE uses the ciphertexts Enc_pk (A) and Enc_pk (b) and private random values in order
#             to obtain encryptions of new values called “masked data”
#     2.      These encryptions are then sent to the CSP, which decrypts and runs a given algorithm
#             on the masked data
# 	3.      This computations output is the “masked model”, a vector w~, which is then sent back
# 	        from the CSP to the MLE
# 	4.      MLE computes the output w^* from w~
#
# """
# #Step 1
# #Step 2
# #Step 3
# #Step 4