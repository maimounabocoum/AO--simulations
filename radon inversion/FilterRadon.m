function FILTER = FilterRadon(f, N , filter , Fc)
% Function written by Maïmouna Bocoum - 28/02/2017 as fully inspired by
% iradon function from matlab packages (type help iradon for more details)
%   
% N : size of fourier transform (N = 2^k) used for radon fourier transform
% order = 2^nextpow2(N) ;
% order = 2*max(64,order)
N = N/2 ;

n = 0:N; 
filtImpResp = zeros(1,N+1); % 'filtImpResp' is the bandlimited ramp's impulse response (values for even n are 0)
filtImpResp(1) = 1/4; % Set the DC term 
filtImpResp(2:2:end) = -1./((pi*n(2:2:end)).^2); % Set the values for odd n
% symmetry of filter around 0 axis :
filtImpResp = [filtImpResp filtImpResp(end-1:-1:2)]; 
filt = 2*real(fft(filtImpResp)); 
%FILTER = filt(1:N+1);

FILTER = abs(f).^0.9 ;

switch filter
    case 'ram-lak'
        % Do nothing
    case 'shepp-logan'
        % be careful not to divide by 0:
        FILTER = FILTER.* (sin(2*pi*f/(2*Fc))./((2*pi*f)/(4*Fc)));  
    case 'cosine'
        %FILTER(2:end) = FILTER(2:end) .* cos(2*pi*f(2:end)/(0.5*2*Fc));
        FILTER = FILTER .* cos(2*pi*f/(4*Fc));
    case 'hamming'
        FILTER = FILTER.* (.54 + .46 * cos(2*pi*f/Fc));
    case 'hann'
        FILTER = FILTER.*(1+cos(2*pi*f./(4*Fc))) / 2;
    otherwise
        error(message('images:iradon:invalidFilter'))
end

%FILTER(isnan(FILTER)) = 0 ;

 FILTER(f > Fc ) = 0;                         % Crop the frequency response
% FILTER = [FILTER' ; FILTER(end-1:-1:2)'];    % Symmetry of the filter

end

