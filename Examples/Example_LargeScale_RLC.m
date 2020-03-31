% -------------------------------------------------------------------------
% This example demonstrates the method for Example 7.3 in the paper
%
% Benner, P., Goyal, P., and Van Dooren, P., Identification of
% Port-Hamiltonian systems from Frequency Response Data. arXiv:1911.00080,
% 2019. 
%
% The model in this example is taken from the paper "S. Gugercin, 
% A. C. Antoulas, A survey of balancing methods for model
% reduction, in: Proc. European Control Conf. ECC 2003, Cambridge, UK,
% 2003".
% -------------------------------------------------------------------------
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Affero General Public License for more details.
%
% You should have received a copy of the GNU Affero General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% Copyright (C) 2019-2020 Peter Benner, Pawan Goyal, Paul Van Dooren
% -------------------------------------------------------------------------

clearvars; clc; close all
addpath(genpath('../Functions/'))

%% Define a system
load data/rlc_serkan200.mat 
E       = eye(size(A)); 
[n,m]   = size(B);

% Define the original transfer function 
H_orig  = @(s) C*((s*E-A)\B) + D;


%% Generate the frequency response data by taking Np in a defined frequency range
Np      = 2*n;
w = logspace(-1,3,Np);
F = zeros(1,Np);
for j = 1:Np
   F(j) = H_orig(1i*w(j));  
end
%% Identification of pH model
tol = 1e-8;
[IdenpHModel, out] = Construct_pHModel(w,F,D,tol);
SingVals_LoewPencil = out.SingVals;

%% Compare the transfer functions of the original and identified systems
w       = logspace(-2,4,500);
sH_orig = zeros(1,length(w));
sH_Loew = zeros(1,length(w));
sH_Iden = zeros(1,length(w));
for j = 1:length(w)
    sH_orig(j) = H_orig(w(j)*1i);
    sH_Loew(j) = out.LoewModel.TF(w(j)*1i);
    sH_Iden(j) = IdenpHModel.TF(w(j)*1i);
end

figure(1)
subplot(2,1,1)
loglog(w,abs(sH_orig),'b',w,abs(sH_Loew),'--g',w,abs(sH_Iden),'-.r')
legend('original','Loewner','Identified pH Model','Location','Best')
subplot(2,1,2)
loglog(w,abs(sH_orig-sH_Loew),'--g',w,abs(sH_orig-sH_Iden),'-.r')
legend('Loewner','Identified pH Model','Location','Best')

%% Plot the Spectral zeros of the original and identified pH system 
% Compute Spectral zero of the original system
GenEign.A     = [zeros(n,n), A, B; A', zeros(n,n) C'; B', C, D + D'];
GenEign.E     = [zeros(n,n) E zeros(n,m); -E, zeros(n,n), zeros(n,m); zeros(m,n), zeros(m,n),zeros(m,m)];

[eig_Val_SZ] = eig(GenEign.A,GenEign.E);
SpecZeros    = diag(eig_Val_SZ); SpecZeros = SpecZeros(1:end-m); % Select finite eigenvalues

pH_SpecZeros = out.SpecZeros; % Spectral zero of the identified pH system

figure(2)
plot(real(SpecZeros), imag(SpecZeros),'r+')
hold all
plot(real(pH_SpecZeros), imag(pH_SpecZeros),'go')
title('Spectral zero')
legend('Using Original Model','Loewner Model','Location','Best')
xlabel('Real part')
ylabel('Imag part')

%% Plot the decay of the singular values of the Loewner pencil
figure(3)
semilogy(SingVals_LoewPencil, 'LineWidth',2)