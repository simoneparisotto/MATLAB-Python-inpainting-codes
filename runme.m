%% MATLAB Codes for the Image Inpainting Problem
%
% Authors:
% Simone Parisotto          (email: sp751 at cam dot ac dot uk)
% Carola-Bibiane Schoenlieb (email: cbs31 at cam dot ac dot uk)
%      
% Address:
% Cambridge Image Analysis
% Centre for Mathematical Sciences
% Wilberforce Road
% CB3 0WA, Cambridge, United Kingdom
%  
% Date:
% September, 2016
%
% Licence: BSD-3-Clause (https://opensource.org/licenses/BSD-3-Clause)
%
clear
close all
clc

addpath ./lib
addpath ./dataset

%% AMLE (Absolute Minimizing Lipschitz Extension) Inpainting
clear
close all

% create the corrupted image with the mask
cleanfilename = 'amle_clean.png';
maskfilename  = 'amle_mask.png';
[u,mask]      = create_image_and_mask(cleanfilename,maskfilename);
imwrite(u,'./dataset/amle_input.png')

% parameters
lambda        = 10^2; 
tol           = 1e-8;
maxiter       = 40000;
dt            = 0.01;

% inpainting
tic
inpainting_amle(u,mask,lambda,tol,maxiter,dt);
toc

%% Harmonic Inpainting
clear
close all

% create the corrupted image with the mask
cleanfilename = 'harmonic_clean.png';
maskfilename  = 'harmonic_mask.png';
[u,mask]      = create_image_and_mask(cleanfilename,maskfilename);
imwrite(u,'./dataset/harmonic_input.png')

% parameters
lambda        = 10;
tol           = 1e-5;
maxiter       = 500;
dt            = 0.1;

% inpainting
tic
inpainting_harmonic(u,mask,lambda,tol,maxiter,dt);
toc

%% Mumford-Shah Inpainting
clear
close all

% create the corrupted image with the mask
cleanfilename = 'mumford_shah_clean.png';
maskfilename  = 'mumford_shah_mask.png';
[u,mask]      = create_image_and_mask(cleanfilename,maskfilename);
imwrite(u,'./dataset/mumford_shah_input.png')

% parameters
maxiter       = 20; 
tol           = 1e-14;
param.lambda  = 10^9;   % weight on data fidelity (should usually be large).
param.alpha   = 1;      % regularisation parameters \alpha.
param.gamma   = 0.5;    % regularisation parameters \gamma.
param.epsilon = 0.05;   % accuracy of Ambrosio-Tortorelli approximation of the edge set.

% inpainting
tic
inpainting_mumford_shah(u,mask,maxiter,tol,param);
toc

%% Cahn-Hilliard Inpainting
clear
close all

% create the corrupted image with the mask
cleanfilename = 'cahn_hilliard_clean.png';
maskfilename  = 'cahn_hilliard_mask.png';
[u,mask]      = create_image_and_mask(cleanfilename,maskfilename);
imwrite(u,'./dataset/cahn_hilliard_input.png')

% parameters
maxiter       = 4000; 
param.epsilon = [100 1];
param.lambda  = 10;
param.dt      = 1;

% inpainting
tic
inpainting_cahn_hilliard(u,mask,maxiter,param);
toc

%% Transport Inpainting
clear
close all

% create the corrupted image with the mask
cleanfilename = 'transport_clean.png';
maskfilename  = 'transport_mask.png';
[u,mask]      = create_image_and_mask(cleanfilename,maskfilename);
imwrite(u,'./dataset/transport_input.png')

% parameters
tol           = 1e-5;
maxiter       = 50;
dt            = 0.1;
param.M       = 40; % number of steps of the inpainting procedure;
param.N       = 2;  % number of steps of the anisotropic diffusion;
param.eps     = 1e-10;

% inpainting
tic
inpainting_transport(u,mask,maxiter,tol,dt,param);
toc
