import numpy as np

branchid = 4
A = np.genfromtxt('%d.csv' % branchid, delimiter=',')
index = A[:,0] % 50 == 0
B = A[index,:]
np.savetxt("%d_bis.csv"%branchid, B, delimiter=",")

