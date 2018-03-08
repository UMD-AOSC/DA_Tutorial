import numpy as np

R = [[1, 2, 3],
     [4, 5, 6],
     [7, 8, 9]]

i = [0,2]
j = i

R = np.array(R)
R2 = R[i,:]
R2 = R2[:,j]

print('R = ')
print(R)
print('R2 = ')
print(R2)
