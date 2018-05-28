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

function inpainting_harmonic(imagefilename,maskfilename,lambda,tol,maxiter,dt)

%% ----------------- CREATE A log FILE WHERE TO STORE RESULTS IN txt FORMAT
logfilename = 'log_harmonic.log';
if exist(logfilename,'file')
    delete(logfilename);
end
fileID = fopen(logfilename,'w');

%% ------------------------------------ IMPORT THE CLEAN INPUT AND THE MASK
iminfo = imfinfo(imagefilename);
input  = im2double(imread(imagefilename));
% check if grayscale/truecolor dimension of image grey/colour
colors = size(input,3);

mask   = im2double(imread(maskfilename));
mask = double(mat2gray(mask)>0);
if size(mask,3)==1 && colors>1
    mask = repmat(mask,[1,1,colors]);
end

%% ---------------------------------------------------- INITIALIZATION OF u
u_start = (~mask).*input + mask;
u_end   = zeros(iminfo.Height,iminfo.Width,colors);

%% ---------------------------------------------- GRID INTERVAL FOR AXIS ij
h1 = 1; h2 = 1; 

%% -------------------------------------------------------------- LAPLACIAN
d2i = toeplitz(sparse([1,1],[1,2],[-2,1]/h1^2,1,iminfo.Height));
d2j = toeplitz(sparse([1,1],[1,2],[-2,1]/h2^2,1,iminfo.Width));
% NEUMANN BOUNDARY CONDITIONS
d2i(1,[1 2])         = [-1 1]/h1;
d2i(end,[end-1 end]) = [1 -1]/h1;
d2j(1,[1 2])         = [-1 1]/h2;
d2j(end,[end-1 end]) = [1 -1]/h2;
% 2D domain LAPLACIAN
L = kron(speye(iminfo.Width),d2i)+kron(d2j,speye(iminfo.Height));

%% ------------------------------------------------------------ FREE MEMORY
clear d2i d2j

%% -------------------------------------------------------------- ALGORITHM
% INITIALIZATION
u            = reshape(u_start,iminfo.Height*iminfo.Width,colors);
f            = reshape(u_start,iminfo.Height*iminfo.Width,colors);
channel_mask = reshape(mask,iminfo.Height*iminfo.Width,colors);

% FOR EACH COLOR CHANNEL
for k=1:colors
    for iter = 1:maxiter
        % COMPUTE NEW SOLUTION
        unew = u(:,k) + dt*(L*u(:,k) + lambda*(~channel_mask(:,k)).*(f(:,k)-u(:,k)));
        
        % COMPUTE EXIT CONDITION
        diff = norm(unew-u(:,k))/norm(unew);
        
        % UPDATE
        u(:,k) = unew;
        
        % TEST EXIT CONDITION
        if diff<tol
            break
        end
    end
    
    % WRITE ON log FILE
    fprintf(fileID,'Channel %d: Iterations: %d, Normalised difference of u: %2.4e\n',k,iter,diff);
    u_end(:,:,k) = reshape(u(:,k),iminfo.Height,iminfo.Width);
end

fclose(fileID);

%%  ---------------------------------------------------- WRITE IMAGE OUTPUT
imwrite(u_start,'masked_harmonic.png')
imwrite(u_end,'output_harmonic.png')

return