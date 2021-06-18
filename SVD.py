import numpy as np 
import math

def comple(arr, valid):
    assert np.shape(arr)[0] == np.shape(arr)[1], "Invalid matrix!"
    iter = 5
    n = arr.shape[0]
    for i in range(valid, n):
        x =  np.random.randn(n,1)
        for j in range(5):
            for k in range(i):
                tmp = arr[:, k].reshape((-1, 1))
                x -= tmp * tmp.T.dot(x)
            x /= np.linalg.norm(x)

        arr[:, i] = x[:, 0]
    return arr            


# Diagonize a real symmetric matrix(n*n), return eigenvalues vector(1*n) and orthogonal matrix(n*n)
def Diag(a):                                        
    row,col = np.shape(a)
    U = np.zeros_like(a,dtype = 'float64')
    S = np.zeros(row,dtype = 'float64')
    iter = 100
    iszero = 0
    for i in range(row):                        # Use the iteration x_{n+1} = ||Ax_{n}|| , x = x / ||x|| to find x with respect to eigenvalue with largest norm 
        x =  np.random.randn(row,1)

        for j in range(iter):                   # Gram-Smith normalization.
            for k in range(i):
                x -= np.reshape( U[:,k] * U[:,k].transpose().dot(x) , (row,1) )

            x = np.dot(a,x)
            if np.linalg.norm(x) < 1e-12:       # The remaining eigenvalues are zeros, use complement operation
                iszero = 1
                comple (U, i)
                break
            else:
                x /= np.linalg.norm(x)

        if not iszero:
            S[i] = x.transpose().dot(a).dot(x)
            U[:,i] = x[:,0]
        else:
            break 
    return S, U

# The svd implementation with A(m*n)
def SVD(A):          
    is_tranposed = False
    if A.shape[0] > A.shape[1]:
        is_tranposed = True
        A = A.T   

    row,col = np.shape(A)
    MU = A.dot(A.transpose())           # Seperately calculate left/right eigenvector matrix

    s,U = Diag(MU)
    # print("----", U.shape)

    S = np.zeros_like(A, dtype = 'float64')  # S is the sigma matrix(m*n) with singular values diagonally. 
    for i in range(row):
        S[i,i] = math.sqrt(s[i])

    MV = U.T.dot(A)
    V = np.zeros( (col, col))
    for i in range(col):
        if i < row and S[i, i] > 1e-12:
            V[:, i] = MV[i,:] / S[i, i]
        else:
            V = comple (V, i)
            break
    if not is_tranposed:
        return U, S, V
    else:
        return V, S.T, U

if __name__ == '__main__':
    a = np.array([[0,0.5,0.25],[0.5,0.25,0.125],[0.25,0.125,0.5*0.125]]) 
    [u,s,v] = SVD(a)
    print(u)
    print(s)
    print(v)
    # for i in range(100):
        # m = np.random.randint(1, 1001)
        # n = np.random.randint(1, 11)
        # a = np.random.rand(m, n) 
        # print(a)
        # print(i + 1)
        # assert np.linalg.norm( u.T.dot(u) - np.eye (a.shape[0])) / np.size(a) < 1e-3, "U is not a valid normal matrix"
        # assert np.linalg.norm( v.T.dot(v) - np.eye (a.shape[1]) ) / np.size(a) < 1e-3, "V is not a valid normal matrix"
        # assert np.linalg.norm( u.dot(s).dot(v.T) - a ) / np.size(a) < 1e-3, "SVD decomposition is incorrect"