function [IdenpHModel, out] = Construct_pHModel(w,F,D, tol)
% -------------------------------------------------------------------------
% This function identifies an underlying port-Hamiltonian realization using frequency
% response data. The methodology is proposed in the paper
%
% Benner, P., Goyal, P., and Van Dooren, P., Identification of
% Port-Hamiltonian systems from Frequency Response Data. arXiv:1911.00080,
% 2019. 
%
% -------------------------------------------------------------------------
% The function inputs are:
% w     --  Frequencies on the j-omega axis, omega > 0,
% F     --  Transfer function values at the frequencies w,
% D     --  Direct feed-through term
% tol (otional)   --  tolerance for the Loewner pencil that determines the order of
%           the identified state-space port-Hamiltonian model. Its default
%           value is 1e-8.
% -------------------------------------------------------------------------
%
% The function outputs are:
% IdenpHModel -- Identified port-Hamiltonian the system matrices and transfer
%                   function
% Out         -- It stores additional information:
%                * An intermediate state-space Loewner model and its
%                  transfer function 
%                * Spectral zeros of the identified port-Hamiltonian
%                  system.
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
% Contact: Pawan Goyal, goyalp@mpi-magdeburg.mpg.de
% -------------------------------------------------------------------------

out = [];
if  nargin < 4
    tol = 1e-8;
end
if nargin < 3
    error('Error! The number of inputs should be 3. Please provide frequencies transfer function values at these frequencies, and the direct feed-through term.')
end
%
Np = length(w);
m  = size(D,1);
[L,sL,~,~,V,W,TL,TR]  = loewner(1i*w,F);

% Matching at infinity and transformation into real realization.
L   = real(TL*L*TR);             
sL  = real(TL*sL*TR - TL*ones(Np,1)*D*ones(1,Np)*TR);     
W    = real(W*TR - ones(1,Np)*D*TR);
V    = real(TL*V - TL*D*ones(Np,1));

% Compression of the realization and construction of minimal realization
[Y,svL, ~]   = svd([L sL]);
[~,~,X]      = svd([L; sL]);
svL          = diag(svL);
r            = sum((svL./svL(1))>tol);
Yr           = Y(:,1:r);
Xr           = X(:,1:r);

% State space model;
El = -Yr'*L*Xr;
Al = -Yr'*sL*Xr;
Bl = Yr'*V;
Cl = W*Xr;

out.LoewModel.E     = El; 
out.LoewModel.A     = Al; 
out.LoewModel.B     = Bl; 
out.LoewModel.C     = Cl;
out.LoewModel.D     = D;
out.LoewModel.TF    = @(s) Cl*((s*El-Al)\Bl) + D;

% Spectral zero of the Loewner Model
%% Compute Spectral zero and directions
GenEign.A     = [zeros(r,r), Al, Bl; Al', zeros(r,r) Cl'; Bl', Cl, D + D'];
GenEign.E     = [zeros(r,r) El zeros(r,m); -El', zeros(r,r), zeros(r,m); zeros(m,r), zeros(m,r),zeros(m,m)];

[eig_Vec_SZLoew, ...
    eig_Val_SZLoew] = eig(GenEign.A,GenEign.E);
SpecZeros            = diag(eig_Val_SZLoew); SpecZeros = SpecZeros(1:end-m); % Select finite eigenvalues
SpecZeros_pos        = SpecZeros(real(SpecZeros)>0); % Select positive eigenvalues

% Extract zero directions: these are the part of the eigenvectors
dir     = zeros(m,2*r+m);
for i = 1:2*r+m
    dir(:,i) = eig_Vec_SZLoew(end-m+1:end,i)/norm(eig_Vec_SZLoew(end-m+1:end,i));
end
dir_pos = dir(:,real(SpecZeros)>0); % select the directions corresponding to positive zero spectral.

%% Now, using the Loewner model, we construct new points that lead us to a pH realization.
data    = zeros(m,length(SpecZeros_pos));
for i = 1:length(SpecZeros_pos)
    data(:,i) = out.LoewModel.TF(SpecZeros_pos(i))*dir_pos(:,i);
end
%% 
% Solve for Loewner and shifted Loewner matrices using corresponding Slyvester form (However, note that there exist analytical solutions for 
% Loewmer and shifted Loewner. But to avoid confusion with notation, we
% take the Sylvester route. Since they are in a reduced-dimension,
% computation cost is not a big issue. 

Lambda  = diag(SpecZeros_pos);
R       = dir_pos;
W       = data;

Loew_RHS    = R'*W + W'*R;
sLoew_RHS   = R'*W*Lambda - Lambda'*W'*R;

% Compute Loewner matrix
Loew = Loew_RHS./((SpecZeros_pos + SpecZeros_pos').');

% Compute Shifted Loewner matrix
sLoew = sLoew_RHS./((SpecZeros_pos + SpecZeros_pos').');

B_iden = -W';
C_iden = -W;

% Shift the shifted-Loewner, B, and C matrices in order to match the
% function at infinity
IdenModel.D     = D;
sLoew_Mod       = sLoew - R'*IdenModel.D*R;
B_iden_Mod      = B_iden - R'*IdenModel.D;
C_iden_Mod      = C_iden + IdenModel.D *R;

% Transfer the realization into real
t           = (1/sqrt(2))*[1 1i; 1 -1i];
T = [];
i = 1;
while i <= length(SpecZeros_pos)
    if abs(imag(SpecZeros_pos(i))) >= 1e-10
        T = blkdiag(T,t);
        i = i + 2;
    else 
        T = blkdiag(T,1);
        i = i + 1;
    end
end

IdenpHModel.E   = real(T'*Loew*T);
IdenpHModel.A   = real(T'*sLoew_Mod*T);
IdenpHModel.B   = real(T'*B_iden_Mod);
IdenpHModel.C   = real(C_iden_Mod*T);
IdenpHModel.D   = D;
IdenpHModel.TF  = @(s) IdenpHModel.C*((s*IdenpHModel.E-IdenpHModel.A)\IdenpHModel.B) + IdenpHModel.D; 

out.SpecZeros = SpecZeros;
out.SingVals  = svL;