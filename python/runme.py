# PYTHON Codes for the Image Inpainting Problem
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

import sys
sys.path.append("./lib")

from ttictoc import tic,toc

import matplotlib.image as mpimg
import inpainting

### AMLE Inpainting

# create the corrupted image with the mask
cleanfilename = './dataset/amle_clean.png';
maskfilename  = './dataset/amle_mask.png';
input,mask    = inpainting.create_image_and_mask(cleanfilename,maskfilename);
mpimg.imsave("./dataset/amle_input.png", input[:,:,0], cmap="gray")

# parameters
fidelity      = 10^2; 
tol           = 1e-8;
maxiter       = 40000;
dt            = 0.01;

# inpainting
tic()
u = inpainting.amle(input,mask,fidelity,tol,maxiter,dt)
print('Elasped time is:',toc())

### Harmonic Inpainting

# create the corrupted image with the mask
cleanfilename = './dataset/harmonic_clean.png';
maskfilename  = './dataset/harmonic_mask.png';
input,mask    = inpainting.create_image_and_mask(cleanfilename,maskfilename);
mpimg.imsave("./dataset/harmonic_input.png", input)

# parameters
fidelity      = 10;
tol           = 1e-5;
maxiter       = 500;
dt            = 0.1;

# inpainting
tic()
u = inpainting.harmonic(input,mask,fidelity,tol,maxiter,dt)
print('Elasped time is:',toc())

### Mumford-Shah Inpainting

# create the corrupted image with the mask
cleanfilename = './dataset/mumford_shah_clean.png';
maskfilename  = './dataset/mumford_shah_mask.png';
input,mask    = inpainting.create_image_and_mask(cleanfilename,maskfilename);
mpimg.imsave("./dataset/mumford_shah_input.png", input)

# parameters
maxiter   = 20; 
tol       = 1e-14;
fidelity  = 10^9;   # weight on data fidelity (should usually be large).
alpha     = 1;      # regularisation parameters \alpha.
gamma     = 0.5;    # regularisation parameters \gamma.
epsilon   = 0.05;   # accuracy of Ambrosio-Tortorelli approximation of the edge set.

# inpainting
tic()
u = inpainting.mumford_shah(input,mask,maxiter,tol,fidelity,alpha,gamma,epsilon);
print('Elasped time is:',toc())

### Cahn-Hilliard Inpainting

# create the corrupted image with the mask
cleanfilename = './dataset/cahn_hilliard_clean.png';
maskfilename  = './dataset/cahn_hilliard_mask.png';
input,mask    = inpainting.create_image_and_mask(cleanfilename,maskfilename);
mpimg.imsave("./dataset/cahn_hilliard_input.png", input)

# parameters
maxiter  = 4000
fidelity = 10
dt       = 1

# inpainting
tic()
u = inpainting.cahn_hilliard(input,mask,maxiter,fidelity,dt)
print('Elasped time is:',toc())

### Transport Inpainting

# create the corrupted image with the mask
cleanfilename = './dataset/transport_clean.png';
maskfilename  = './dataset/transport_mask.png';
input,mask    = inpainting.create_image_and_mask(cleanfilename,maskfilename);
mpimg.imsave("./dataset/transport_input.png", input)

# parameters
maxiter          = 50;
tol              = 1e-5;
dt               = 0.1;
iter_inpainting  = 40; # number of steps of the inpainting procedure;
iter_anisotropic = 2;  # number of steps of the anisotropic diffusion;
epsilon          = 1e-10;

# inpainting
tic()
u = inpainting.transport(input,mask,maxiter,tol,dt,iter_inpainting,iter_anisotropic,epsilon)
print('Elasped time is:',toc())
