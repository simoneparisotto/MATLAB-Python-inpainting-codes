%% AMLE Inpainting 
% part of "MATLAB Codes for the Image Inpainting Problem"
%
% Usage:
% u = inpainting_amle(u,mask,lambda,tol,maxiter,dt)
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

function u = inpainting_amle(input,mask,lambda,tol,maxiter,dt)

[M,N,C] = size(input);

%% GRADIENT
% GRID INTERVAL FOR AXIS ij
h1 = 1; h2 = 1;

% average upper (u_{i+1,j} - u_{ij})/hx  and lower (u_{ij}-u_{i-1,j})/hx
d1i_forward  = spdiags([-ones(M,1),ones(M,1)],[0,1],M,M)/h1;
d1i_backward = spdiags([-ones(M,1),ones(M,1)],[-1,0],M,M)/h1;
% average upper (u_{i,j+1} - u_{ij})/hy  and lower (u_{ij}-u_{i,j-1})/hy
d1j_forward  = spdiags([-ones(N,1),ones(N,1)],[0,1],N,N)/h2;
d1j_backward = spdiags([-ones(N,1),ones(N,1)],[-1,0],N,N)/h2;

% BACKWARD WITHOUT BOUNDARY CONDITIONS (FORWARD NOT OF INTEREST HERE)
% DD1i_forward  = kron(speye(N),d1i_forward);
DD1i_backward = kron(speye(N),d1i_backward);
% DD1j_forward  = kron(d1j_forward,speye(M));
DD1j_backward = kron(d1j_backward,speye(M));

% NEUMANN BOUNDARY CONDITION
d1i_forward(end,:) = 0; d1i_backward(1,:)  = 0;
d1j_forward(end,:) = 0; d1j_backward(1,:)  = 0;

% FORWARD AND BACKWARD WITH BOUNDARY CONDITIONS
D1i_forward  = kron(speye(N),d1i_forward);
D1i_backward = kron(speye(N),d1i_backward);
D1j_forward  = kron(d1j_forward,speye(M));
D1j_backward = kron(d1j_backward,speye(M));

% CENTERED WITH BOUNDARY CONDITIONS
D1i_centered = (D1i_forward+D1i_backward)/2;
D1j_centered = (D1j_forward+D1j_backward)/2;

%% FREE MEMORY
clear d1i_forward d1i_backward d1j_forward d1j_backward

%% INPAINTING ALGORITHM
% INITIALIZATION
u     = reshape(input,M*N,C);
input = reshape(input,M*N,C);
mask  = reshape(mask,M*N,C);

v     = zeros(M*N,2);

% FOR EACH COLOR CHANNEL
for c = 1:C

    % ITERATION
    for iter = 1:maxiter
        ux = D1i_forward*u(:,c); % forward differences along i
        uy = D1j_forward*u(:,c); % forward differences along j
    
        % second derivatives
        uxx = DD1i_backward*ux;
        uxy = DD1j_backward*ux;
        uyx = DD1i_backward*uy;
        uyy = DD1j_backward*uy;
    
        % create direction field Du/|Du| with central differences
        v(:,1) = D1i_centered*u(:,c);
        v(:,2) = D1j_centered*u(:,c);
        % normalize the direction field
        v = bsxfun(@rdivide,v,sqrt(sum(v.^2,2)));
        v(isnan(v)) = 0;
    
        % CORE ITERATION
        unew = u(:,c) + dt*(uxx.*v(:,1).^2+uyy.*v(:,2).^2 + (uxy+uyx) .* (v(:,1).*v(:,2)) + lambda*mask(:,c).*(input(:,c)-u(:,c)));
    
        % COMPUTE EXIT CONDITION
        diff = norm(unew-u(:,c))/norm(unew);
    
        % UPDATE
        u(:,c) = unew;
    
        % TEST EXIT CONDITION
        if diff<tol
            break
        end
    end
    
end

%% GET THE 2D SOLUTION
u = reshape(u,M,N,C);

%% WRITE IMAGE OUTPUT
imwrite(u,'./results/amle_output.png')

return