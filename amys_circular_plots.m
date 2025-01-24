%Read the table & assign varible names

path='/Users/amgi2571/Desktop/personal_workspace/xyz_aspect/liv_24.csv';

array=readtable(path);
data=array.Var1;

%convert data (degrees) into radians
rad_data=circ_ang2rad(data);

%plot data in circular histogram https://www.mathworks.com/help/matlab/ref/polarhistogram.html#d126e1260832
plot=polarhistogram(rad_data,20,'FaceColor','#004d40','FaceAlpha',.8,'LineWidth',1.3);
title('Amys Polar Plot');


%calculate mean resultant vector length (this test or Rayleigh test??)
vectorlength=circ_r(rad_data)


%Ominibus calculation
[p_value, m]=circ_otest(rad_data)






%%%%%FUNCTIONS%%%%%%
%degrees to radians
function alpha = circ_ang2rad(alpha)
alpha = alpha * pi /180;
end

%mean resultant vector length
function r = circ_r(alpha, w, d, dim)
% r = circ_r(alpha, w, d, dim)

if nargin < 4
  dim = find(size(alpha) > 1, 1, 'first');
  if isempty(dim)
    dim = 1;
  end
end

if nargin < 2 || isempty(w) 
  % if no specific weighting has been specified
  % assume no binning has taken place
	w = ones(size(alpha));
else
  if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1) 
    error('Input dimensions do not match');
  end 
end

if nargin < 3 || isempty(d)
  % per default do not apply correct for binned data
  d = 0;
end

% compute weighted sum of cos and sin of angles
r = sum(w.*exp(1i*alpha),dim);

% obtain length 
r = abs(r)./sum(w,dim);

% for data with known spacing, apply correction factor to correct for bias
% in the estimation of r (see Zar, p. 601, equ. 26.16)
if d ~= 0
  c = d/2/sin(d/2);
  r = c*r;
end
end

%omnibus calculation
function [pval, m] = circ_otest(alpha, sz, w)

if size(alpha,2) > size(alpha,1)
	alpha = alpha';
   
end

if nargin < 2 || isempty(sz)
  sz = circ_ang2rad(1);
  
end

if nargin < 3
  w = ones(size(alpha));
else
  if length(alpha)~=length(w)
    error('Input length does not match.')
  end
  w =w(:);  
end

alpha = mod(alpha,2*pi);
n = sum(w);
dg = 0:sz:pi;

m1 = zeros(size(dg));
m2 = zeros(size(dg));
for i=1:length(dg)
  m1(i) = sum((alpha > dg(i) & alpha < pi + dg(i)).*w);    
  m2(i) = n - m1(i);
end
m = min(min([m1;m2]));

if n > 50
  % approximation by Ajne (1968)
  A = pi*sqrt(n) / 2 / (n-2*m);
  pval = sqrt(2*pi) / A * exp(-pi^2/8/A^2);
else
  % exact formula by Hodges (1955)
  % pval = 2^(1-n) * (n-2*m) * nchoosek(n,m);  % revised below for numerical stability
  pval = exp((1-n)*log(2) + log(n-2*m) + gammaln(n+1) - gammaln(m+1) - gammaln(n-m+1));
end
end
