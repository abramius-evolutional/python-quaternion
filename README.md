## Example

	import numpy as np
	
	q = Quaternion(hi=np.pi/2, vector=[1., 0., 0.])
	v = np.array([0., 1., 1.])
	
	# Rotate vector
	print q.rotateVector(v)
	# Rotate basis for vector
	print q.rotateBasis(v)
	# Getting transition matrix
	print q.getMatrix()
