# operators module
# part of "PYTHON Codes for the Image Inpainting Problem"
#
# Authors:
# Simone Parisotto          (email: sp751 at cam dot ac dot uk)
# Carola-Bibiane Schoenlieb (email: cbs31 at cam dot ac dot uk)
#      
# Address:
# Cambridge Image Analysis
# Centre for Mathematical Sciences
# Wilberforce Road
# CB3 0WA, Cambridge, United Kingdom
#  
# Date:
# October, 2018
#
# Licence: BSD-3-Clause (https://opensource.org/licenses/BSD-3-Clause)
#

    
import scipy
from scipy.sparse import spdiags
from scipy.linalg import toeplitz

def GRAD_ij(X,hi,hj,direction,boundary):
    """
    Compute the GRADIENT OPERATOR of an image along i and j dimension
    """
    M, N = X.shape       
        
    if direction=='backward':
       D1 = spdiags([scipy.ones(M),-scipy.ones(M)],[-1,0],M,M)/hi
       D2 = spdiags([scipy.ones(N),-scipy.ones(N)],[-1,0],N,N)/hj
    elif direction=='central':
       D1 = spdiags([-scipy.ones(M),scipy.ones(M)],[-1,1],M,M)/(2*hi)
       D2 = spdiags([-scipy.ones(N),scipy.ones(N)],[-1,1],N,N)/(2*hj)
    elif direction=='forward':
       D1 = spdiags([-scipy.ones(M),scipy.ones(M)],[0,1],M,M)/hi
       D2 = spdiags([-scipy.ones(N),scipy.ones(N)],[0,1],N,N)/hj

    D1 = D1.tolil()
    D2 = D2.tolil()
    
    if boundary=='neumann':
        if direction=='backward':
            D1[0,:] = 0
            D2[0,:] = 0
        elif direction=='central':
            D1[(0,M-1),:] = 0
            D2[(0,M-1),:] = 0
        elif direction=='forward':
            D1[M-1,:] = 0
            D2[N-1,:] = 0   
    elif boundary=='periodic':
        if direction=='forward':
            D1[0,1]     = 0
            D2[0,1]     = 0
            D1[M-1,M-2] = 0
            D2[N-1,N-2] = 0
            D1[0,M-1]   = -1
            D2[0,N-1]   = -1
            D1[M-1,0]   = 1
            D2[N-1,0]   = 1
        if direction=='backward':
            D1[0,1]     = 0
            D2[0,1]     = 0
            D1[M-1,M-2] = 0
            D2[N-1,N-2] = 0
            D1[0,M-1]   = 1
            D2[0,N-1]   = 1
            D1[M-1,0]   = -1
            D2[N-1,0]   = -1
            
    #D1 = scipy.sparse.kron( scipy.sparse.identity(N),D1 )
    #D2 = scipy.sparse.kron( D2,scipy.sparse.identity(M) )
    D1 = scipy.sparse.kron( D1,scipy.sparse.identity(N) )
    D2 = scipy.sparse.kron( scipy.sparse.identity(M),D2 )
    
    return (D1,D2)

def DIV_ij(D1,D2):
    """
    Compute the DIVERGENCE OPERATOR of a vector field along i and j dimension
    """
    
    D1T = -D1.T
    D2T = -D2.T
    
    return (D1T,D2T)

def LAP_ij(X,hi,hj,direction,boundary):
    """
    Compute the LAPLACIAN OPERATOR of an image
    """
    D1,D2   = GRAD_ij(X,hi,hj,direction,boundary)
    D1T,D2T = DIV_ij(D1,D2)
    L       = D1T.dot(D1) + D2T.dot(D2)
    
    
    return L

def laplacian_ij(X,hi,hj):
    """
    Compute the LAPLACIAN OPERATOR of an image
    """
    M,N = X.shape 
    
    D2i = (toeplitz(scipy.append(scipy.array([-2,1]),scipy.zeros((M-2,1))))/(hi**2));
    D2j = (toeplitz(scipy.append(scipy.array([-2,1]),scipy.zeros((N-2,1))))/(hj**2));

    D2i[0,1]     = 2/(hi**2)
    D2i[M-2,M-1] = 2/(hi**2)
    D2j[0,1]     = 2/(hj**2)
    D2j[N-2,N-1] = 2/(hj**2)
    
    #L = scipy.sparse.kron(scipy.sparse.identity(N),D2i) + scipy.sparse.kron(D2j,scipy.sparse.identity(M))
    L = scipy.sparse.kron(D2i,scipy.sparse.identity(N)) + scipy.sparse.kron(scipy.sparse.identity(M),D2j)
    return L
