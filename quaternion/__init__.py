import numpy as np

class Quaternion:
    def _vectorProduct(self, v1, v2):
        '''vectors format: numpy.matrix([ [0.0],[0.0],[0.0] ]) or numpy.array([ [0.0],[0.0],[0.0] ])'''
        result = np.matrix( [ [0.0],[0.0],[0.0] ] )
        result[0,0] = v1[1,0] * v2[2,0] - v1[2,0] * v2[1,0]
        result[1,0] = v1[2,0] * v2[0,0] - v1[0,0] * v2[2,0]
        result[2,0] = v1[0,0] * v2[1,0] - v1[1,0] * v2[0,0]
        return result

    def _mixedProduct(self, v1, v2, v3):
        mart = np.matrix( [ [[0.0],[0.0],[0.0]], [[0.0],[0.0],[0.0]], [[0.0],[0.0],[0.0]] ] )
        mart[0,0], mart[0,1], mart[0,2] = v1[0,0], v1[1,0], v1[2,0]
        mart[1,0], mart[1,1], mart[1,2] = v2[0,0], v2[1,0], v2[2,0]
        mart[2,0], mart[2,1], mart[2,2] = v3[0,0], v3[1,0], v3[2,0]
        return np.linalg.det(mart)
    def _absVector(self, vector):
        '''vector format: numpy.matrix([ [0.0],[0.0],[0.0] ]) or numpy.array([ [0.0],[0.0],[0.0] ])'''
        return ( ( vector[0,0] )**2 + ( vector[1,0] )**2 + ( vector[2,0] )**2 ) ** 0.5

    def _prodQuaternions(self, q1, q2):
        Q = Quaternion()
        Q.scalar = q1.scalar * q2.scalar - np.vdot(np.array(q1.vector),np.array(q2.vector))
        Q.vector = q1.scalar*q2.vector + q2.scalar*q1.vector + self._vectorProduct(q1.vector, q2.vector)
        return Q

    def _conjugate(self, q):
        Q = Quaternion()
        Q.scalar = q.scalar
        Q.vector = q.vector * (-1.)
        return Q

    def _conjugate(self):
        Q = Quaternion()
        Q.scalar = self.scalar
        Q.vector = self.vector * (-1.)
        return Q

    def _norma(self, q):
        return q.scalar**2 + np.vdot(q.vector,q.vector)

    def _norma(self):
        return self.scalar**2 + np.vdot(self.vector,self.vector)

    def _inv(self, q):
        Q = Quaternion()
        _norma = self._norma(q)
        sopr = self._conjugate(q)
        if _norma != 0:
            Q.scalar = sopr.scalar / _norma
            Q.vector = sopr.vector / _norma
        else:
            raise 'error: zero _norma Quaternion'
        return Q

    def _inv(self):
        Q = Quaternion()
        _norma = self._norma()
        sopr = self._conjugate()
        if _norma != 0:
            Q.scalar = sopr.scalar / _norma
            Q.vector = sopr.vector / _norma
        else:
            raise 'error: zero _norma Quaternion'
        return Q

    def rotateVector(self, vector):
        '''np.array([0, 0, 0])'''
        Q = Quaternion()
        Q.scalar = 0.
        Q.vector = np.matrix([[vector[0]], [vector[1]], [vector[2]]])
        sopr = self._conjugate()
        result = self._prodQuaternions( self, self._prodQuaternions( Q, sopr ) ).vector
        return np.array([result[0,0], result[1, 0], result[2, 0]])

    def rotateBasis(self, vector):
        '''np.array([0, 0, 0])'''
        Q = Quaternion()
        Q.scalar = 0.
        Q.vector = np.matrix([[vector[0]], [vector[1]], [vector[2]]])
        sopr = self._conjugate()
        result = self._prodQuaternions( sopr, self._prodQuaternions( Q, self ) ).vector
        return np.array([result[0,0], result[1, 0], result[2, 0]])

    def getMatrix(self):
        M = np.matrix([ [0.0,0,0], [0,0,0], [0,0,0] ])
        L0 = self.scalar
        L1 = self.vector[0,0]
        L2 = self.vector[1,0]
        L3 = self.vector[2,0]
        M[0,0], M[0,1], M[0,2] = L0**2+L1**2-L2**2-L3**2, 2*(L0*L3+L1*L2), 2*(-L0*L2+L1*L3)
        M[1,0], M[1,1], M[1,2] = 2*(-L0*L3+L2*L1), L0**2+L2**2-L3**2-L1**2, 2*(L0*L1+L2*L3)
        M[2,0], M[2,1], M[2,2] = 2*(L0*L2+L3*L1), 2*(-L0*L1+L3*L2), L0**2+L3**2-L1**2-L2**2
        return M

    def __init__(self, hi=np.pi, vector=[0.0, 0.0, 0.0]):
        '''hi - vector of the final rotation, e - direction of the axis of rotation
           vector format: np.array([0, 0, 0])'''
        self.vector = np.matrix([[0.0],[0.0],[0.0]])
        e = np.matrix([[vector[0]], [vector[1]], [vector[2]]])
        if self._absVector(e) != 0:
            ee = e / self._absVector(e)
        else:
            ee = e
        self.scalar = np.cos(hi/2)
        self.vector = ee * np.sin(hi/2)
