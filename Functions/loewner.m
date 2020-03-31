function [LL,sLL,mu,la,V,W,TL,TR] = loewner(s,F,varargin)
% LOEWNER   Loewner and shifted-Loewner matrices associated with a data set.
%
%  [LL,SLL] = LOEWNER(S,F) constructs the Loewner matrix LL and the
%  shifted-Loewner matrix SLL associated with the values F at points S,
%  where S is an N-vector, and F is either an N-vector for the scalar case,
%  or an P-by-M-by-N array for the case when F represents a collection of N
%  P-by-M matrices.
%
%  For an even number of points S, N = 2*K, LL and SLL are K-by-K matrices.
%  For an odd number of points S, N = 2*K+1, they are (K+1)-by-K matrices.
%
%  [LL,SLL] = LOEWNER(S,F,'real') ensures that LL and SLL have only real
%  entries. For a complex (matrix) value F(k), its complex conjugate pair 
%  (conj(S(k)),conj(F(k))) is also added to the data set (S,F),
%  increasing its dimension accordingly.
%  Example: S = [1; 1j];      F = [2; -2j] becomes
%           S = [1; 1j; -1j]; F = [2; -2j; 2j].
%
%  LL and sLL are obtained by first partitioning the data (S,F) into
%  the LEFT subset (MU,V) and RIGHT subset (LA,W)
%
%     MU(1), LA(1), MU(2), LA(2), ... = S(1), S(2), S(3), S(4), ... 
%      V(1),  W(1),  V(2),  W(2), ... = F(1), F(2), F(3), F(4), ... .
%
%  For the SCALAR case, when F is a vector, use the definitions
%
%      LL(i,j) = (       V(i)-W(j)       )/( MU(i)-LA(j) )
%     SLL(i,j) = ( MU(i)*V(i)-W(j)*LA(j) )/( MU(i)-LA(j) ).
%  
%  For the MATRIX case, when F is a P-by-M-by-N array, partition F as
%  before into two subsets of left P-by-M matrices V(:,:,i) and right 
%  P-by-M matrices W(:,:,j). Then, compute the LEFT tangential 1-by-M
%  vectors V(i,:) using random tangential directions L(i,:)
%
%     L(i,1:P) = randn(1,P);     V(i,1:M) = L(i,1:P)*V(1:P,1:M,i)
%
%  and the RIGHT tangential P-by-1 vectors W(:,j) using random tangential
%  directions R(:,j)
%
%     R(1:M,j) = randn(M,1);     W(1:P,j) = W(1:P,1:M,j)*R(1:M,j),
%
%  and use the tangential definitions
%
%   LL(i,j) = (       V(i,:)*R(:,j)-L(i,:)*W(:,j)       )/( MU(i)-LA(j) )
%  SLL(i,j) = ( MU(i)*V(i,:)*R(:,j)-L(i,:)*W(:,j)*LA(j) )/( MU(i)-LA(j) ).
%
%  [LL,SLL,MU,LA,V,W] = LOEWNER(...) returns the quantities associated
%  with the partitioned S and F. MU and LA are column vectors, V is a tall
%  matrix with M columns, W is a wide matrix with P rows. 
%
%  References:
%  	A. J. Mayo, A. C. Antoulas, A framework for the solution of the 
%     generalized realization problem, Linear Algebra and its Applications
%     425, pp. 634?662 (2007). (DOI: 10.1016/j.laa.2007.03.008)
%
%     S. Lefteriu, A. C. Antoulas, A New Approach to Modeling Multiport
%     Systems From Frequency-Domain Data, IEEE Transactions on 
%     Computer-Aided Design of Integrated Circuits and Systems,
%     Vol. 29, No. 1, pp. 14-27 (2010). (DOI: 10.1109/TCAD.2009.2034500)

%  DISCLAIMER:
%  This software package was written for teaching and testing purposes, and
%  is distributed "AS IS", without any written or implied warranties.
%  
%  LICENSE:
%  This software is free. You can modify and/or redistribute it under the
%  terms of the GNU General Public License www.gnu.org/licenses/gpl.html .
%
%  Antonio Cosmin Ionita
%  Rice University
%  2013/05/06

% check if the function inputs are valid
if nargin < 2 || nargin > 3
   error('loewner:InconsistentData','Wrong number of input arguments.');
end
if any([~isnumeric(s), ~isnumeric(F(:)), isempty(s), isempty(F), ...
      any(~isfinite(s)), any(~isfinite(F(:)))])
   error('loewner:InconsistentData','Data set is not numeric.');
end
if ~isvector(s)
   error('loewner:InconsistentData:s','Points are not in a vector.');
end

% ensure the points s are in a column-vector
s = s(:);
N = size(s,1);

% check the dimension of the function values F
sizeF = size(F);
if length(sizeF) == 2
   F = F(:);
   P = 1;            M = 1;            NF = max(sizeF); 
elseif length(sizeF) == 3
   P = sizeF(1);     M = sizeF(2);     NF = sizeF(3); 
else
   error('loewner:InconsistentData:F','Values are not in a vector/array.');
end

% check if number of points and values are the same
if N ~= NF
   error('loewner:InconsistentData', ...
         'The number of points S and values F is not the same.');
end

% check for the 'real' option
rflag = 1;
% if nargin == 3
%    if strcmp(varargin{1},'real')
%       rflag = 1;
%    else
%       error('loewner:InconsistentData', ...
%       'Use loewner(s,F,''real'') to get matrices with real entries.');
%    end
% end

%
% Partition the points and values.
%
idx1 = 1:2:N;                          % mu1, mu2, ... = s1, s3, ...
idx2 = 2:2:N;                          % la1, la2, ... = s2, s4, ...
mu = s(idx1);                          % left points
la = s(idx2);                          % right points
if P == 1 && M == 1                    % case of SCALAR values
   V = F(idx1);                        % left values
   W = F(idx2).';                      % right value
else                                   % case of MATRIX values
   % need tangential directions to handle multiple inputs and outpus
   L = randn(length(mu),P);            % left tangential directions
   R = randn(M,length(la));            % right tangential directions
%   L = ones(length(mu),P);            % left tangential directions
%   R = ones(M,length(la));            % right tangential directions  
   for i = 1:length(mu)
      V(i,:) = L(i,:)*F(:,:,idx1(i));  % left tangential values
   end
   for j = 1:length(la)
      W(:,j) = F(:,:,idx2(j))*R(:,j);	% right tangential values
   end
end

% If Loewner matrices with real entries are requested, append the
% complex conjugate values to the left and right subsets.

if rflag == 1
	% LEFT subset: separate the real entries
   idxRe = (imag(mu) == 0);
   muRe = mu(idxRe);
   if sum(idxRe) ~= length(mu)
      % append the complex entries and their complex conjugates
      idxIm = (imag(mu) ~= 0);
      muIm(1:2:2*sum(idxIm),1) = mu(idxIm);
      muIm(2:2:2*sum(idxIm),1) = conj(mu(idxIm));
      VIm(1:2:2*sum(idxIm),:) = V(idxIm,:);
      VIm(2:2:2*sum(idxIm),:) = conj(V(idxIm,:));
      mu = [muRe; muIm];
      V = [V(idxRe,:); VIm];
      if P > 1 || M > 1
         LIm(1:2:2*sum(idxIm),:) = L(idxIm,:);
         LIm(2:2:2*sum(idxIm),:) = conj(L(idxIm,:));
         L = [L(idxRe,:); LIm];
      end
   end
	% RIGHT subset: separate the real entries
   idxRe = (imag(la) == 0);
   laRe = la(idxRe);
   if sum(idxRe) ~= length(la)
      % append the complex entries and their complex conjugates
      idxIm = (imag(la) ~= 0);
      laIm(1:2:2*sum(idxIm),1) = la(idxIm);
      laIm(2:2:2*sum(idxIm),1) = conj(la(idxIm));
      WIm(:,1:2:2*sum(idxIm)) = W(:,idxIm);
      WIm(:,2:2:2*sum(idxIm)) = conj(W(:,idxIm));
      la = [laRe; laIm];
      W = [W(:,idxRe), WIm];
      if P > 1 || M > 1
         RIm(:,1:2:2*sum(idxIm)) = R(:,idxIm);
         RIm(:,2:2:2*sum(idxIm)) = conj(R(:,idxIm));
         R = [R(:,idxRe), RIm];
      end
   end
end

%
% Form the Loewner and shifted-Loewner matrices.
%

% apply the definitions
eNL = ones(length(mu),1);
eNR = ones(length(la),1);

den = mu*eNR.'-eNL*la.';                           % the denominators

if P == 1 && M == 1                                % scalar case
   LL = ( V*eNR.'-eNL*W ) ./ den;                  % Loewner
   sLL = ( (mu.*V)*eNR.'-eNL*(W.*la.') ) ./ den;   % shifted-Loewner
else                                               % matrix case
   LL = ( V*R-L*W ) ./ den;
   sLL = ( diag(mu)*V*R-L*W*diag(la) ) ./ den;
end

%
% apply coordinate transformation to obtain matrices with real entries
%
% if rflag == 1
   TL = eye(length(muRe));
   TR = eye(length(laRe));   
   for i = 1:length(muIm)/2,  TL = blkdiag(TL,1/sqrt(2)*[1 1;-1j 1j]);  end
   for i = 1:length(laIm)/2,  TR = blkdiag(TR,1/sqrt(2)*[1 1j;1 -1j]);  end
% return
%    transform and discard possible small imaginary parts (due to roundoff)
%    LL = real(TL*LL*TR);             % discard possible imaginary parts
%    sLL = real(TL*sLL*TR);           %  in the order of machine precission
%    the left and right subsets in the new coordinate system
%    W = real(W*TR);
%    V = real(TL*V);
%   Dla = TR'*diag(la)*TR;
%   Dmu = TL*diag(mu)*TL';
end
