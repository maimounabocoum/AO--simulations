% define erbium absorption and emission cross sections
% 
% taking the very scientific approach of visual inspection from a plot
function [abint,emint] = ErSpectrum(wlint)

if nargin<1,
    wlint = 1450:0.5:1640; 
end
if wlint(1)<1
    wlint=wlint*1e9;
end

pk = 0.6e-24;
wl = [980 1450:5:1600, 1610:10:1640];
ab = (pk/100)* ... 
    [28/0.6, ...                                        % 980
        14, 20, 25, 30, 35, 40, 44, 47, 48, 50, 52, ... % 1450 to 1500
        53, 56, 61, 72, 93, 100, 78, 60, 54, 45, ...    % 1505 vto 1550
        40, 34, 26, 19, 14, 10, 8, 7, 5, 4, ...         % 1555 to 1600
        3, 2, 1, 0];                                    % 1610 to 1640
%     
%         ...
%     [0.28, 0.025, 0.05, 0.15, 0.27, 0.24, 0.20, 0.21... % 1460 to 1520
%     0.25, 0.4, 0.6, 0.39, 0.29, 0.24, ... % 1525 to 1550
%     0.14, 0.06, 0.04, 0.03, 0.02, 0.015, 0.01, 0.002, 0.001]; % 1560 to 1640
em = (pk/100)* ... 
    [0/0.6, ...                                         % 980
        4, 4, 5, 7, 10, 12, 15, 17, 20, 23, 27, ...     % 1450 to 1500
        30, 34, 43, 55, 80, 95, 82, 71, 69, 66, ...     % 1505 to 1550
        64, 60, 52, 42, 33, 27, 23, 21, 20, 18, ...     % 1555 to 1600
        14, 8, 0, 0];                                   % 1610 to 1640
% 
% em = (pk/100)*[0, 0.001, 0.02, 0.04, 0.1, 0.12, 0.14, 0.17... % 1460 to 1520
%     0.2, 0.4, 0.6, 0.42, 0.35, 0.3, ... % 1525 to 1550
%     0.24, 0.13, 0.09, 0.1, 0.09, 0.08, 0.07, 0.04, 0.02]; % 1560 to 1640

abint = interp1(wl,ab,wlint,'cubic');
emint = interp1(wl,em,wlint,'cubic');
% figure;plot(wlint,abint,'-b',wlint,emint,':r'); grid on; legend('absorption','emission');

