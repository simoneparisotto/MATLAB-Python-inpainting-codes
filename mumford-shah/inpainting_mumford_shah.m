%% MATLAB Codes for the Image Inpainting Problem
%  Copyright (c) 2016, Simone Parisotto and Carola-Bibiane Schoenlieb
%  All rights reserved.
% 
%  Redistribution and use in source and binary forms, with or without 
%  modification, are permitted provided that the following conditions are met:
% 
%  1. Redistributions of source code must retain the above copyright notice,
%     this list of conditions and the following disclaimer.
% 
%  2. Redistributions in binary form must reproduce the above copyright 
%     notice, this list of conditions and the following disclaimer in the 
%     documentation and/or other materials provided with the distribution.
% 
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
%  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
%  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
%  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
%  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
%  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED 
%  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
%  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
%  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
%  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
%  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%  Authors:
%  Simone Parisotto (email: sp751 at cam dot ac dot uk)
%  Carola-Bibiane Schoenlieb (email: cbs31 at cam dot ac dot uk)
%      
%  Address:
%  Cambridge Image Analysis
%  Centre for Mathematical Sciences
%  Wilberforce Road
%  Cambridge CB3 0WA
%  United Kingdom
%  
%  Date:
%  September, 2016
%%

function inpainting_mumford_shah(imagefilename,maskfilename,maxiter,tol,param)
% Inpainting with the Mumford-Shah image model and Ambrosio-Tortorelli.
% For a given grey value image ustart with image domain \Omega and
% inpainting  domain (damaged part) D we want to reconstruct an image u from f
% by solving
% %%%%%%%%%%%%%% MINIMISATION PROBLEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% u = argmin_u  (\frac{\alpha}{2} \int_\Omega \chi^2 |\nabla u|^2 dx   %
%               + \beta \int_\Omega \left(\epsilon |\nabla \chi|^2     %
%               + \frac{(1-\chi)^2}{4\epsilon} \right) dx              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The above minimisation problem is solved iteratively via alternating
% solutions of the Euler-Lagrange equations for u and \chi.

%% ----------------- CREATE A log FILE WHERE TO STORE RESULTS IN txt FORMAT
logfilename = 'log_mumford_shah.log';
if exist(logfilename,'file')
    delete(logfilename);
end
fileID = fopen(logfilename,'w');

%% ------------------------------------ IMPORT THE CLEAN INPUT AND THE MASK
iminfo = imfinfo(imagefilename);
input  = im2double(imread(imagefilename));
% check if grayscale/truecolor dimension of image grey/colour
colors = size(input,3);

mask = im2double(imread(maskfilename));
mask = double(mat2gray(mask)==0); % indicator function of the intact image
if size(mask,3)==1 && colors>1
    mask = repmat(mask,[1,1,colors]);
end

%% ---------------------------------------------- GRID INTERVAL FOR AXIS ij
h1 = 1/(size(input,1)+1); h2 = 1/(size(input,2)+1);
N = iminfo.Height*iminfo.Width; % number of pixels

%% --------------------------------------------------------------- GRADIENT
% FORWARD AND BACKWARD
d1i_forward  = spdiags([-ones(iminfo.Height,1),ones(iminfo.Height,1)],[0,1],iminfo.Height,iminfo.Height)/h1;
d1j_forward  = spdiags([-ones(iminfo.Width,1),ones(iminfo.Width,1)],[0,1],iminfo.Width,iminfo.Width)/h2;
d1i_backward = spdiags([-ones(iminfo.Height,1),ones(iminfo.Height,1)],[-1,0],iminfo.Height,iminfo.Height)/h1;
d1j_backward = spdiags([-ones(iminfo.Width,1),ones(iminfo.Width,1)],[-1,0],iminfo.Width,iminfo.Width)/h2;
% PERIODIC BOUNDARY CONDITIONS
d1i_forward(end,:) = 0;
d1j_forward(end,:) = 0;
%d1i_forward(end,[1 end]) = [1 -1]/h1;
%d1j_forward(end,[1 end]) = [1 -1]/h2;
d1i_backward(1,:) = 0;
d1j_backward(1,:) = 0;
%d1i_backward(1,[1 end]) = [-1 1]/h1;
%d1j_backward(1,[1 end]) = [-1 1]/h2;

matrices.Dif  = kron(speye(iminfo.Width),d1i_forward);
matrices.Dib  = kron(speye(iminfo.Width),d1i_backward);
matrices.Djf  = kron(d1j_forward,speye(iminfo.Height));
matrices.Djb  = kron(d1j_backward,speye(iminfo.Height));

% CENTRAL
matrices.Dic = (matrices.Dif+matrices.Dib)/2;
matrices.Djc = (matrices.Djf+matrices.Djb)/2;

%% -------------------------------------------------------------- LAPLACIAN
d2i = toeplitz(sparse([1,1],[1,2],[-2,1],1,iminfo.Height))/h1^2;
d2j = toeplitz(sparse([1,1],[1,2],[-2,1],1,iminfo.Width))/h2^2;
% PERIODIC BOUNDARY CONDITIONS
d2i(1,2)       = 2/h1^2;
d2i(end,end-1) = 2/h1^2;
d2j(1,2)       = 2/h2^2;
d2j(end,end-1) = 2/h2^2;
% 2D domain LAPLACIAN
matrices.LAP = (kron(speye(iminfo.Width),d2i)+kron(d2j,speye(iminfo.Height)));

%% ------------------------------------------------------------ FREE MEMORY
clear d1i_forward dji_forward d1i_backward dji_backward
clear d2i d2j

%% -------------------------------------------------------------- ALGORITHM
% (1) INITIALIZATION: u^0 = 0, z^0 = 0, solve for k=1,2,...
mask         = double(mask);
u0           = double(input);
noise        = mat2gray(randn(size(u0)));
u0(~mask)    = noise(~mask);

channel_mask = reshape(mask,[N,colors]);
u            = reshape(u0,[N,colors]);
chi          = (channel_mask);

param.lambda = param.lambda*(channel_mask);

rhsL = (param.lambda/param.gamma).*reshape(u0,[N,colors]);
rhsM = ones(N,colors); 

% FOR EACH COLOR CHANNEL
for k=1:colors
    
    % ITERATION
    for iter = 1:maxiter
  
        
        % SOLVE EULER-LAGRANGE EQUATION FOR \chi
        % i.e M(u^{k-1},\chi^k) = 1.
        % M is a linear operator acting on \chi and reads
        % M(u,.) = 1+\frac{2\epsilon\alpha}{\beta} |\nabla u|^2
        %           - 4\epsilon^2\Delta.
        % Solved via inversion of the linear operators.
        M        = matrixM(param,N,matrices,u(:,k));
        chinew   = M\rhsM(:,k);
        diff_chi = norm(chinew-chi(:,k))/norm(chinew);
        chi(:,k) = chinew;
        clear chinew
        
        % SOLVE EULER-LAGRANGE EQUATION FOR u: 
        % i.e. L(\chi^k,u^k)= \alpha \chi_{\Omega\setminus D} \cdot ustart.
        % L is a linear operator acting on u and reads
        % L(\chi,.) = -\div(\chi_\epsilon^2\nabla)
        %             + \alpha\chi_{\Omega\setminus D}.
        % Solved via inversion of the linear operators.
        L      = matrixL(param,N,matrices,chi(:,k));
        unew   = L\rhsL(:,k);
        diff_u = norm(unew-u(:,k))/norm(unew);
        u(:,k) = unew;
        clear unew
        
        % WRITE ON log FILE
        fprintf(fileID,'Channel %d, normalised difference of u: %2.4e, normalised difference of chi: %2.4e\n',k,diff_u,diff_chi);
        
        % TEST EXIT CONDITION
        if diff_u<tol
            break
        end    
        
    end    
end

fclose(fileID);

%% ---------------------------------------------------- GET THE 2D SOLUTION
u_end   = reshape(u,[iminfo.Height,iminfo.Width,colors]);
chi_end = mat2gray(reshape(chi,[iminfo.Height,iminfo.Width,colors]));

%% ---------------------------------------------------- WRITE IMAGE OUTPUTS
imwrite(u0,'masked_mumford_shah.png')
imwrite(u_end,'output_mumford_shah.png')
imwrite(chi_end,'levels_mumford_shah.png')

return

%% ---------------------------------------------------- AUXILIARY FUNCTIONS
function M = matrixM(param,N,matrices,u)
% Definition of (\nabla u)^2:
nablau2 = (matrices.Dic*u).^2 + (matrices.Djc*u).^2;

M = speye(N)...
    + 2 * param.epsilon * param.gamma/param.alpha * spdiags(nablau2,0,N,N)...
    - (4*param.epsilon^2)*matrices.LAP;
return

function L = matrixL(param,N,matrices,chi)
% Definition of the nonlinear diffusion weighted by \chi^2:

z  = (chi).^2+param.epsilon^2; % coefficient of nonlinear diffusion

zx = matrices.Dic*z;
zy = matrices.Djc*z;
Z  = spdiags(z,0,N,N);
Zx = spdiags(zx,0,N,N);
Zy = spdiags(zy,0,N,N);

NonlinearDelta = Z*matrices.LAP + Zx*matrices.Dic + Zy*matrices.Djc;

L =  -NonlinearDelta + spdiags(param.lambda/param.gamma,0,N,N);

return
