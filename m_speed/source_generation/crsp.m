%**************************************************************************
% Copyright (C) 2012 The SPEED FOUNDATION
% Author: SPEED-Team (POLIMI)
%         Politecnico di Milano 
%         P.zza Leonardo da Vinci, 32 
%         20133 Milano 
%         Italy                                        
%
% This file is part of SPEED.
%
% SPEED is free software; you can redistribute it and/or modify it
% under the terms of the GNU Affero General Public License as
% published by the Free Software Foundation, either version 3 of the
% License, or (at your option) any later version.
%
% SPEED is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Affero General Public License for more details.
%
% You should have received a copy of the GNU Affero General Public License
% along with SPEED.  If not, see <http://www.gnu.org/licenses/>.
%**************************************************************************

function [crsp] = crsp(L,W,ax,ay,dx,dy,slipi,mu,cov,eta,type)

% INPUT PARAMETERS
% L = fault length (along-strike) in m
% W = fault width (down-dip) in m
% ax = along-strike correlation length (in m) 
% ay = down dip correlation length  (in m)
% slipi: matrix of SLIP DISTRIBUTION (obtained after interpolation from 
% Martin Mai databse : http://www.seismo.ethz.ch/srcmod/) 
% dx, dy = spatial sampling steps of the slip interpolation grid 
% mu = mean value of the variabile to be randomly perturbated 
% cov = covariance coeff. 
% eta = correlation coefficient of the random variable with the slip
% distribution matrix
% type = type of Power Spectral Density PSD 

% OUTPUT 
% crsp = matrix correlated random source variable (dimension: ndipxnstr)


%disp('not yet implemented')
%stop


% input data for source parameter correlated random perturbations
sigma = cov.*mu; % standard deviation

% slip distribution (format as in SCRMOD database, http://www.seismo.ethz.ch/srcmod/) 
dim_scrmod = size(slipi);
M = dim_scrmod(2);
N = dim_scrmod(1);
% normalization of slip distribution matrix
% it is assumed that S is a uniform distribution 
S = (slipi)./abs(max(max(slipi))); 

% wave number discretization 
dkx = 1/(M*dx); 
dky = 1/(N*dy); 
kx_Ny = 1/2/dx; 
ky_Ny = 1/2/dy; 
kx = [0:dkx:(M-1)*dkx]; 
ky = [0:dky:(N-1)*dky];

%==========================================================================%
% 2D Power Spectral Density-PSD ( see Frankel & Clayton, 1986 JGR
%                                 see Mai & Beroza, 2002, JGR     ) 
%==========================================================================%
for k=1:M/2+1
   for l=1:N/2+1
       if type == 0  % GAUSSIAN (GS) CORRELATION FUNCTION 
           PSD(l,k) = ax*ay/2*exp(-1/4*(kx(k)^2*ax^2+ky(l)^2*ay^2)); 
       elseif type == 1  % EXPONENTIAL (EX) CORRELATION FUNCTION 
           PSD(l,k) = ax*ay/((1+kx(k)^2*ax^2+ky(l)^2*ay^2)^1.5); 
       elseif type == 2        % SELF-SIMILAR VON KARMAN (VK) CORRELATION FUNCTION     
           PSD(l,k) = ax*ay/(1+kx(k)^2*ax^2+ky(l)^2*ay^2);      
       end
   end
end

%=========================================================================%
% 2D Gaussian (normal) distributed white noise matrice 
%=========================================================================%
xi1 = zeros(N,M); xi2 = zeros(N,M); xi = zeros(N,M); 
xi1 = random('Normal',0,1,N,M);
xi2 = eta.*S + sqrt(1-eta^2).*xi1;
xi = xi2; 
%FFT of random distribution xi 
XI = fft2(xi,N,M);


%=========================================================================%
% SPECTRAL FILTERING METHOD (SPECTRAL REPRESENTATION ALGORITHM)  
%=========================================================================%
for k = 1:M/2+1
   for l = 1:N/2+1
       if (k~=1 || l~=1) 
             F(l,k) = sqrt(PSD(l,k))/sqrt(2).*(real(XI(l,k))+sqrt(-1)*imag(XI(l,k)));
       else
           F(l,k) = 0.0; 
       end
       
       if k==1 
          k0 = 1; 
       else 
          k0 = M-k+2; 
       end
       if l==1 
          l0 = 1; 
       else 
          l0 = N-l+2; 
       end        
       F(l0,k0) = conj(F(l,k)); 
   end
end
F(1,M/2+1) = real(sqrt(PSD(1,M/2+1))/sqrt(2).*(real(XI(1,M/2+1))+sqrt(-1)*imag(XI(1,M/2+1)))); 
F(N/2+1,1) = real(sqrt(PSD(N/2+1,1))/sqrt(2).*(real(XI(N/2+1,1))+sqrt(-1)*imag(XI(N/2+1,1))));
F(N/2+1,M/2+1) = real(sqrt(PSD(N/2+1,M/2+1))/sqrt(2).*(real(XI(N/2+1,M/2+1))+sqrt(-1)*imag(XI(N/2+1,M/2+1))));
for k=2:M/2
   for l=2:N/2 
       F(N-l+2,k) = conj(F(l,k)); 
       F(l,M-k+2) = conj(F(l,k));
   end
end

fxy = ifft2(F,N,M);
fxy = real(fxy);
fxy = fxy./abs(max(max(fxy)));
crsp = fxy*sigma + mu; 


return 


 


