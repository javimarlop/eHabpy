import scipy.ndimage as nd
import numpy as np

a = np.array([[0,0,1,1,0,0], # this is the similar area map
                  [0,0,0,1,0,0],
                  [1,1,0,0,1,0],
                  [0,0,0,1,0,0]])

print a

s = nd.generate_binary_structure(2,2)

print s

labeled_array, num_features = nd.label(a, structure=s)

print(num_features)

print(labeled_array)

sum1 = sum(labeled_array[labeled_array==1]/1)
sum2 = sum(labeled_array[labeled_array==2]/2)

print sum1
print sum2

coms = nd.measurements.center_of_mass(a,labeled_array,[1,2])

print coms


