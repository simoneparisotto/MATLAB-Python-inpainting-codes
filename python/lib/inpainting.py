# inpainting module
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


from __future__ import division

import numpy as np
import scipy
from scipy.sparse import spdiags
from scipy.sparse import linalg
import pyfftw
import matplotlib.image as mpimg
import cv2 as cv
from math import pi
import sys
sys.path.append("./lib")
import operators as opy


### Create image and mask
def create_image_and_mask(imagefilename,maskfilename):

    # import a clean input to be corrupted with the mask 
    input = mpimg.imread(imagefilename)
    
    if input.ndim==3:
        M,N,C = input.shape
    else:
        M,N = input.shape
        C = 1
    
    # import the mask of the inpainting domain
    # mask = 1 intact part
    # mask = 0 missing domain
    mask  = scipy.float64((mpimg.imread(maskfilename) == 1));
    
    if (input.ndim==3) & (mask.ndim<3):
        mask = np.repeat(mask[:, :, np.newaxis], C, axis=2)
    
    if C==1:
        input = scipy.expand_dims(input, axis=2)
        mask  = scipy.expand_dims(mask, axis=2)
        
    # create the image with the missin domain:
    noise = scipy.rand(M,N,C)      
    u     = mask*input + (1-mask)*noise;
    
    return (u,mask)

def create_kernel_derivatives():
    
    # KERNELS
    # forward i,j
    Kfi = scipy.zeros((3,3))
    Kfi[1,1] = -1
    Kfi[2,1] = 1
    Kfj = scipy.zeros((3,3))
    Kfj[1,1] = -1
    Kfj[1,2] = 1

    # backward i,j
    Kbi = scipy.zeros((3,3))
    Kbi[1,1] = 1
    Kbi[0,1] = -1
    Kbj = scipy.zeros((3,3))
    Kbj[1,1] = 1
    Kbj[1,0] = -1

    # centred i,j
    Kci = scipy.zeros((3,3))
    Kci[0,1] = -0.5
    Kci[2,1] = 0.5
    Kcj = scipy.zeros((3,3))
    Kcj[1,2] = 0.5
    Kcj[1,0] = -0.5
    
    return (Kfi,Kfj,Kbi,Kbj,Kci,Kcj)

### AMLE Inpainting
def amle(input,mask,fidelity,tol,maxiter,dt):
    
    if input.ndim==3:
        M,N,C = input.shape
    else:
        M,N = input.shape
        C = 1 
    
    # KERNELS of derivatives    
    Kfi,Kfj,Kbi,Kbj,Kci,Kcj = create_kernel_derivatives()
    
    # INITIALISATION

    u = input.copy()
    v = scipy.zeros((M,N,2));

    # ITERATION
    for c in range(0,C):
        
        for iter in range(0,maxiter):
        
            ux = cv.filter2D(u[:,:,c],-1, Kfi) # forward differences along i
            uy = cv.filter2D(u[:,:,c],-1, Kfj) # forward differences along j
        
            # second derivatives
            uxx = cv.filter2D(ux,-1, Kbi)
            uxy = cv.filter2D(ux,-1, Kbj)
            uyx = cv.filter2D(uy,-1, Kbi)
            uyy = cv.filter2D(uy,-1, Kbj)
    
            # create direction field Du/|Du| with central differences
            v[:,:,0] = cv.filter2D(u[:,:,c],-1, Kci)
            v[:,:,1] = cv.filter2D(u[:,:,c],-1, Kcj)
        
            # normalize the direction field
            dennormal = scipy.sqrt( scipy.sum(v**2,axis=2) + 1e-15)
            v[:,:,0] = v[:,:,0]/dennormal
            v[:,:,1] = v[:,:,1]/dennormal
    
            # CORE ITERATION
            unew = u[:,:,c] + dt*( uxx*v[:,:,0]**2 + uyy*v[:,:,1]**2 + (uxy+uyx)*(v[:,:,0]*v[:,:,1]) + fidelity*mask[:,:,c]*( input[:,:,c]-u[:,:,c] ) );

            # exit condition
            diff_u = np.linalg.norm(unew.reshape(M*N,1)-u[:,:,c].reshape(M*N,1),2)/np.linalg.norm(unew.reshape(M*N,1),2); 

            # update
            u[:,:,c] = unew
     
            # test exit condition
            if diff_u<tol:
                break

    if C==1:
        mpimg.imsave("./results/amle_output.png", u[:,:,0],cmap="gray")
    elif C==3:
        mpimg.imsave("./results/amle_output.png", u)
        
    return u

### Harmonic Inpainting
def harmonic(input,mask,fidelity,tol,maxiter,dt):

    if input.ndim==3:
        M,N,C = input.shape
    else:
        M,N = input.shape
        C = 1

    u = input.copy()

    for c in range(0,C):
    
        for iter in range(0,maxiter):
    
            # COMPUTE NEW SOLUTION
            laplacian = cv.Laplacian(u[:,:,c],cv.CV_64F)
            unew = u[:,:,c] + dt*( laplacian + fidelity * mask[:,:,c] * (input[:,:,c]-u[:,:,c]) )
    
            # exit condition
            diff_u = np.linalg.norm(unew.reshape(M*N,1)-u[:,:,c].reshape(M*N,1),2)/np.linalg.norm(unew.reshape(M*N,1),2); 

            # update
            u[:,:,c] = unew
     
            # test exit condition
            if diff_u<tol:
                break
    
    
    mpimg.imsave("./results/harmonic_output.png", u)
    
    return u


### Mumford-Shah Inpainting

def build_matrixM(u,alpha,gamma,epsilon,Dic,Djc,L):
    
    N, = u.shape
    
    # Definition of \nabla u)^2:
    nablau2 = Dic.dot( u.ravel() )**2 + Djc.dot( u.ravel() )**2
    
    M = scipy.sparse.eye(N,N) + (2*epsilon*gamma/alpha)*scipy.sparse.spdiags(nablau2,0,N,N) - 4*(epsilon*2)*L;
    
    return M

def build_matrixL(chi,FIDELITY,gamma,epsilon,Dic,Djc,L):
    
    N, = chi.shape
    
    # Definition of the nonlinear diffusion weighted by \chi^2:
    z  = chi**2 + epsilon**2; # coefficient of nonlinear diffusion

    zx = Dic.dot(z)
    zy = Djc.dot(z)

    Z  = scipy.sparse.spdiags(z,0,N,N);
    Zx = scipy.sparse.spdiags(zx,0,N,N);
    Zy = scipy.sparse.spdiags(zy,0,N,N);

    NonlinearDelta = Z.dot(L) + Zx.dot(Dic) + Zy.dot(Djc);

    L =  -NonlinearDelta + scipy.sparse.spdiags(FIDELITY/gamma,0,N,N);

    return L    

def mumford_shah(input,mask,maxiter,tol,fidelity,alpha,gamma,epsilon):
    
    if input.ndim==3:
        M,N,C = input.shape
    else:
        M,N = input.shape
        C = 1
    
    u        = (input.copy()).reshape(M*N,C)
    chi      = mask.copy().reshape(M*N,C)
    FIDELITY = fidelity*chi;
    rhsL     = (FIDELITY/gamma)*u
    rhsM     = scipy.ones((M*N,C)); 
    
    hi = 1
    hj = 1
    
    Dfi,Dfj = opy.GRAD_ij(input[:,:,0],hi,hj,'forward','periodic')
    Dbi,Dbj = opy.GRAD_ij(input[:,:,0],hi,hj,'backward','periodic')
    
    Dci     = (Dfi + Dbi)/2
    Dcj     = (Dfj + Dbj)/2

    L = opy.laplacian_ij(input[:,:,0],hi,hj)
   
    # ITERATION
    for c in range(0,C):
    
        # FOR EACH COLOR CHANNEL
        for iter in range(0,maxiter):
      
            # SOLVE EULER-LAGRANGE EQUATION FOR \chi
            # i.e M(u^{c-1},\chi^c) = 1.
            # M is a linear operator acting on \chi and reads
            # M(u,.) = 1+\frac{2\epsilon\alpha}{\beta} |\nabla u|^2
            #           - 4\epsilon^2\Delta.
            # Solved via inversion of the linear operators.
            MM       = build_matrixM(u[:,c],alpha,gamma,epsilon,Dci,Dcj,L)
            chinew   = scipy.sparse.linalg.spsolve(MM,rhsM[:,c]);
            #diff_chi = np.linalg.norm(chinew-chi[:,c],2)/np.linalg.norm(chinew,2);
            chi[:,c] = chinew;
        
            # SOLVE EULER-LAGRANGE EQUATION FOR u: 
            # i.e. L(\chi^c,u^c)= \alpha \chi_{\Omega\setminus D} \cdot ustart.
            # L is a linear operator acting on u and reads
            # L(\chi,.) = -\div(\chi_\epsilon^2\nabla)
            #             + \alpha\chi_{\Omega\setminus D}.
            # Solved via inversion of the linear operators.
            LL     = build_matrixL(chi[:,c],FIDELITY[:,c],gamma,epsilon,Dci,Dcj,L)
            unew   = scipy.sparse.linalg.spsolve(LL,rhsL[:,c]);
            diff_u = np.linalg.norm(unew-u[:,c],2)/np.linalg.norm(unew,2);
            u[:,c] = unew;
            
            # TEST EXIT CONDITION
            if diff_u<tol:
                break

    u   = u.reshape(M,N,C);
    chi = chi.reshape(M,N,C);

    # WRITE IMAGE OUTPUTS
    if C==1:
        mpimg.imsave("./results/mumford_shah_output.png", u[:,:,0],cmap="gray")
    elif C==3:
        mpimg.imsave("./results/mumford_shah_output.png", u)

    # imwrite(chi,'./results/mumford_shah_levels_output.png')

    return u

### Cahn-Hilliard Inpainting
def cahn_hilliard(input,mask,maxiter,fidelity,dt):
    
    if input.ndim==3:
        M,N,C = input.shape
    else:
        M,N = input.shape
        C = 1
    hi = 1
    hj = 1
    
    epsilon  = scipy.array([100,1])

    swap = round(maxiter/2)
    ep1  = scipy.array( epsilon[0]*scipy.ones( (swap,1) ) )
    ep2  = scipy.array( epsilon[1]*scipy.ones( (maxiter-swap,1  ) ) )
    ep   = scipy.concatenate([ep1,ep2])
    c1   = 1.0 / epsilon[1]

    # FIDELITY MASK
    FIDELITY = fidelity * mask[:,:,1]

    # Diagonalize the Laplace Operator by: Lu + uL => D QuQ + QuQ D, where 
    # Q is nonsingular, the matrix of eigenvectors of L and D is a diagonal matrix.
    # We have to compute QuQ. This we can do in a fast way by using the fft-transform:

    l1 = scipy.array( 2*(scipy.cos( 2*scipy.array(range(0,M))*pi / M ) - 1) )
    l2 = scipy.array( 2*(scipy.cos( 2*scipy.array(range(0,N))*pi / N ) - 1) )
    Lambda1 = spdiags(l1,0,M,M)/(hi**2);
    Lambda2 = spdiags(l2,0,N,N)/(hj**2);

    Denominator1 = Lambda1.dot(scipy.ones((M,N))) 
    Denominator2 = Lambda2.T.dot(scipy.ones((N,M))).T
    Denominator  = Denominator1 + Denominator2

    u = scipy.ones((M,N,C))

    for c in range(0,C):

        uu              = input[:,:,c]
        u_hat           = pyfftw.interfaces.numpy_fft.fft2( uu )
        fidelity_u0_hat = pyfftw.interfaces.numpy_fft.fft2( FIDELITY * uu )

        for iter in range(0,maxiter):
        
            fidelity_u_hat = pyfftw.interfaces.numpy_fft.fft2( FIDELITY * uu);
        
            Fprime_hat = pyfftw.interfaces.numpy_fft.fft2( 2*(2*uu**3 - 3*uu**2 + uu) )
        
            # CH-inpainting
            u1   = (1 + dt*FIDELITY) * u_hat
            u2   = -( dt/epsilon[1] ) * Denominator * u_hat
    
            u3   = (dt / ep[iter])*Denominator * Fprime_hat
        
            u4   = dt * ( fidelity_u0_hat-fidelity_u_hat )

            uden = 1.0 + dt*( FIDELITY + ep[iter]*(Denominator**2) - (Denominator / epsilon[1]) )
    
            u_hat =  (u1 + u2 + u3 + u4)/uden;
        
            uu = pyfftw.interfaces.numpy_fft.ifft2(u_hat)
    
        u[:,:,c] = uu.real
        
    if C==1:
        mpimg.imsave("./results/cahn_hilliard_output.png", u[:,:,0],cmap="gray")
    elif C==3:
        mpimg.imsave("./results/cahn_hilliard_output.png", u)

    return u
   
### Transport Inpainting
    
# anisotropic diffusion
def anisodiff(u,dt,eps,geps,iterations,Kfi,Kfj,Kbi,Kbj,Kci,Kcj):

    # M. Bertalmio, "Processing of flat and non-flat image information on 
    # arbitrary manifolds using Partial Differential Equations", PhD Thesis, 2001.
    
    M,N = u.shape
    
    Kii = Kfi - Kbi
    Kjj = Kfj - Kbj
    ux  = cv.filter2D(u,-1,Kci)
    uy  = cv.filter2D(u,-1,Kcj)
    uxx = cv.filter2D(u,-1,Kii)
    uyy = cv.filter2D(u,-1,Kjj)
    uxy = cv.filter2D(cv.filter2D(u,-1,Kci),-1,Kcj)
    
    for i in range(0,iterations): 
        squared_normgrad = ux**2 + uy**2 + eps;
        u = u + dt*geps*(uyy * (ux**2) + uxx * (uy**2) - 2*ux*uy*uxy) / squared_normgrad;

    return u

def transport(input,mask,maxiter,tol,dt,iter_inpainting,iter_anisotropic,epsilon):
        
    if input.ndim==3:
        M,N,C = input.shape
    else:
        M,N = input.shape
        C = 1

    # KERNELS of derivatives    
    Kfi,Kfj,Kbi,Kbj,Kci,Kcj = create_kernel_derivatives()

    # INITIALISATION
    laplacian_gepsilon = scipy.zeros((M,N,C))
    Dlaplacian         = scipy.zeros((M,N,2))
    laplacian          = scipy.zeros((M,N,C))
    update             = scipy.zeros((M,N,C))
    zeros              = scipy.zeros((M,N))
    normal             = scipy.zeros((M,N,2))
    LAMBDAEPS          = scipy.zeros((M,N,C));

    u                  = input.copy()
    un                 = input.copy()
    MASKEPS            = 1-mask.copy()
    gepsilon           = 1-mask.copy();

    # MORPHOLOGIC ELEMENT for the small epsilon-diffusion strip around the inpainting domain:
    SE = cv.getStructuringElement(cv.MORPH_ELLIPSE,(13,13))

    for c in range(0,C):
    
        # INTERPOLATE g_{epsilon} with a few steps of linear diffusion within the strip
        LAMBDAEPS[:,:,c] = cv.dilate(MASKEPS[:,:,c], SE, iterations=1) - MASKEPS[:,:,c]
        for t in range(0,5):
            laplacian_gepsilon[:,:,c] = cv.filter2D(gepsilon[:,:,c],-1, Kfi-Kbi+Kfj-Kbj)
            gepsilon[:,:,c] = gepsilon[:,:,c] + dt*laplacian_gepsilon[:,:,c] + LAMBDAEPS[:,:,c]*( MASKEPS[:,:,c]-gepsilon[:,:,c] );
        
        # ANISOTROPIC DIFFUSION PREPROCESSING STEP
        u[:,:,c] = anisodiff(u[:,:,c],dt,epsilon,gepsilon[:,:,c],1,Kfi,Kfj,Kbi,Kbj,Kci,Kcj);

        # INPAINTING
        for outer_iter in range(0,maxiter):
    
            for iter in range(0,iter_inpainting):
            
                laplacian[:,:,c]   = cv.filter2D(u[:,:,c],-1, Kfi-Kbi+Kfj-Kbj)
                Dlaplacian[:,:,0]  = cv.filter2D(laplacian[:,:,c],-1, Kci)
                Dlaplacian[:,:,1]  = cv.filter2D(laplacian[:,:,c],-1, Kcj)
                normal[:,:,0]      = cv.filter2D(u[:,:,c],-1, Kci)
                normal[:,:,1]      = cv.filter2D(u[:,:,c],-1, Kcj)

                dennormal = scipy.sqrt( scipy.sum(normal**2,axis=2) + epsilon )
        
                normal[:,:,0] = scipy.divide(normal[:,:,0],dennormal)
                normal[:,:,1] = scipy.divide(normal[:,:,1],dennormal)
        
                beta    = Dlaplacian[:,:,0]*(-normal[:,:,1]) + Dlaplacian[:,:,1]*normal[:,:,0]
                betapos = beta>0
        
                uxf = cv.filter2D(u[:,:,c],-1, Kfi)
                uxb = cv.filter2D(u[:,:,c],-1, Kbi)
                uyf = cv.filter2D(u[:,:,c],-1, Kfj)
                uyb = cv.filter2D(u[:,:,c],-1, Kbj)
       
                sl1 = scipy.array([scipy.minimum(uxb,zeros), scipy.maximum(uxf,zeros), scipy.minimum(uyb,zeros), scipy.maximum(uyf,zeros)])
                sl2 = scipy.array([scipy.maximum(uxb,zeros), scipy.minimum(uxf,zeros), scipy.maximum(uyb,zeros), scipy.minimum(uyf,zeros)])
                slopelim1 = scipy.sum(sl1**2,axis=0)
                slopelim2 = scipy.sum(sl2**2,axis=0)
            
                slopelim      = betapos*scipy.sqrt(slopelim1) + (1-betapos)*scipy.sqrt(slopelim2)
                update[:,:,c] = beta * slopelim;
            
                # UPADTE ONLY PIXELS INSIDE THE INPAINTING DOMAIN (BY MASK)
                u[:,:,c] = u[:,:,c] + dt * (1-mask[:,:,c])*update[:,:,c];
        
            # param.N steps of anisotropic diffusion
            un[:,:,c] = anisodiff(u[:,:,c],dt,epsilon,gepsilon[:,:,c],iter_anisotropic,Kfi,Kfj,Kbi,Kbj,Kci,Kcj);

            # exit condition
            diff_u = np.linalg.norm(un[:,:,c].reshape(M*N,1)-u[:,:,c].reshape(M*N,1),2)/np.linalg.norm(un[:,:,c].reshape(M*N,1),2); 

            # update
            u[:,:,c] = (1-mask[:,:,c])*un[:,:,c] + mask[:,:,c]*u[:,:,c];
     
            # test exit condition
            if diff_u<tol:
                break
             
    if C==1:
        mpimg.imsave("./results/cahn_hilliard_output.png", u[:,:,0],cmap="gray")
    elif C==3:
        mpimg.imsave("./results/transport_output.png", u)
          
    return u