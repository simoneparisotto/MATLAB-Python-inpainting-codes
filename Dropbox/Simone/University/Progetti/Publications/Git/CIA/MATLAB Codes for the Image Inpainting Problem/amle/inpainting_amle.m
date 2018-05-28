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

function inpainting_amle(imagefilename,lambda,tol,maxiter,dt)

%% ----------------- CREATE A log FILE WHERE TO STORE RESULTS IN txt FORMAT
logfilename = 'log_amle.log';
if exist(logfilename,'file')
    delete(logfilename);
end
fileID = fopen(logfilename,'w');

%% ------------------------------------ IMPORT THE CLEAN INPUT AND THE MASK
iminfo = imfinfo(imagefilename);
input  = im2double(imread(imagefilename));
mask   = lambda*double(input>0);

%% ---------------------------------------------------- INITIALIZATION OF u
u = rand([iminfo.Height,iminfo.Width]);
u(mask>0) = input(mask>0);

%% ---------------------------------------------- GRID INTERVAL FOR AXIS ij
h1 = 1; h2 = 1;

%% --------------------------------------------------------------- GRADIENT
% average upper (u_{i+1,j} - u_{ij})/hx  and lower (u_{ij}-u_{i-1,j})/hx
d1i_forward  = spdiags([-ones(iminfo.Height,1),ones(iminfo.Height,1)],[0,1],iminfo.Height,iminfo.Height)/h1;
d1i_backward = spdiags([-ones(iminfo.Height,1),ones(iminfo.Height,1)],[-1,0],iminfo.Height,iminfo.Height)/h1;
% average upper (u_{i,j+1} - u_{ij})/hy  and lower (u_{ij}-u_{i,j-1})/hy
d1j_forward  = spdiags([-ones(iminfo.Width,1),ones(iminfo.Width,1)],[0,1],iminfo.Width,iminfo.Width)/h2;
d1j_backward = spdiags([-ones(iminfo.Width,1),ones(iminfo.Width,1)],[-1,0],iminfo.Width,iminfo.Width)/h2;

% BACKWARD WITHOUT BOUNDARY CONDITIONS (FORWARD NOT OF INTEREST HERE)
% DD1i_forward  = kron(speye(iminfo.Width),d1i_forward);
DD1i_backward = kron(speye(iminfo.Width),d1i_backward);
% DD1j_forward  = kron(d1j_forward,speye(iminfo.Height));
DD1j_backward = kron(d1j_backward,speye(iminfo.Height));

% NEUMANN BOUNDARY CONDITION
d1i_forward(end,:) = 0; d1i_backward(1,:)  = 0;
d1j_forward(end,:) = 0; d1j_backward(1,:)  = 0;

% FORWARD AND BACKWARD WITH BOUNDARY CONDITIONS
D1i_forward  = kron(speye(iminfo.Width),d1i_forward);
D1i_backward = kron(speye(iminfo.Width),d1i_backward);
D1j_forward  = kron(d1j_forward,speye(iminfo.Height));
D1j_backward = kron(d1j_backward,speye(iminfo.Height));

% CENTERED WITH BOUNDARY CONDITIONS
D1i_centered = (D1i_forward+D1i_backward)/2;
D1j_centered = (D1j_forward+D1j_backward)/2;

%% -------------------------------------------------------------- LAPLACIAN
d2i = toeplitz(sparse([1,1],[1,2],[-2,1]/h1^2,1,iminfo.Height));
d2j = toeplitz(sparse([1,1],[1,2],[-2,1]/h2^2,1,iminfo.Width));
% NEUMANN BOUNDARY CONDITIONS
d2i(1,[1 2])         = [-1 1]/h1;
d2i(end,[end-1 end]) = [1 -1]/h1;
d2j(1,[1 2])         = [-1 1]/h2;
d2j(end,[end-1 end]) = [1 -1]/h2;
% 2D domain LAPLACIAN
D2i = kron(speye(iminfo.Width),d2i);
D2j = kron(d2j,speye(iminfo.Height));
L = D2i+D2j;

%% ------------------------------------------------------------ FREE MEMORY
clear d1i_forward d1i_backward d1j_forward d1j_backward
clear d2i d2j L

%% -------------------------------------------------------------- ALGORITHM
% INITIALIZATION
u     = u(:);
input = input(:);
mask  = mask(:);
v     = zeros(iminfo.Height*iminfo.Width,2);

% ITERATION
for iter=1:maxiter
    ux = D1i_forward*u; % forward differences along i
    uy = D1j_forward*u; % forward differences along j
    
    % second derivatives
    uxx = DD1i_backward*ux;
    uxy = DD1j_backward*ux;
    uyx = DD1i_backward*uy;
    uyy = DD1j_backward*uy;
    
    % create direction field Du/|Du| with central differences
    v(:,1) = (D1i_centered*u);
    v(:,2) = (D1j_centered*u);
    % normalize the direction field
    v = bsxfun(@rdivide,v,sqrt(sum(v.^2,2)));
    v(isnan(v)) = 0;
    
    % CORE ITERATION
    unew = u + dt*(uxx.*v(:,1).^2+uyy.*v(:,2).^2 + (uxy+uyx) .* (v(:,1).*v(:,2)) + mask.*(input-u));
    
    % COMPUTE EXIT CONDITION
    if ~mod(iter-1,1000)
        diff = norm(unew-u)/norm(unew);
    end
    
    % UPDATE
    u = unew;
    
    % TEST EXIT CONDITION
    if diff<tol
        break
    end
end

% WRITE ON log FILE
fprintf(fileID,'Iterations: %d, Normalised difference of u: %2.4e\n',iter,diff);
fclose(fileID);

%% ---------------------------------------------------- GET THE 2D SOLUTION
u = reshape(u,iminfo.Height,iminfo.Width);

%%  ---------------------------------------------------- WRITE IMAGE OUTPUT
imwrite(u,'output_amle.png')

return