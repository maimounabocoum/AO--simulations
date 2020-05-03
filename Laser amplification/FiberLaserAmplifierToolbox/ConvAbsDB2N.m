function N = ConvAbsDB2N(alpha,sigma)
% convert from absorption in dB to concentration using knowledge of
% absorption coefficients
%
% assume input alpha is given in dB/m as is common in specification sheets
% 
% output will be concentration in parts per m^3

% default values
if nargin<1, alpha = 10; end
if nargin<2, sigma = 0.2e-24; end

% absorbed power using alpha
z = 1;
Pz = 10^(z*alpha/10); % power absorbed from 0 to z

% absorbed power using N and sigma
N = log(Pz) / (sigma*z);

