import pandas as pd
import numpy as np


# D = pd.read_csv('Course_Notes/data/solubility.csv')
# np.save('solubility.npy', np.array(D.Solubility))
D = np.load('Course_Notes/data/solubility.npy')
print(D)
np.savetxt('Course_Notes/data/solubility.txt', D)


# import scipy.io
# # _url = 'http://quantum-machine.org/data/qm7b.mat'
# # mat = scipy.io.loadmat(_url)

# mat = scipy.io.loadmat('Course_Notes/data/qm7b.mat')
# for index, key in enumerate(mat):
#     print(index, key, np.shape(mat[key]))
# print(mat["names"][0])
