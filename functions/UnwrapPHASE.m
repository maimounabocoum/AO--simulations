function [UnwrapANGLE] = UnwrapPHASE(ANGLE,Ix,Iy)
%UNWRAPPHASE Summary of this function goes here
%   Detailed explanation goes here
% Ix: reference point for unwraping
% Iy: center point for unwraping


PhaseRef = ANGLE(Iy,Ix);


UnwrapANGLE = 0*ANGLE;
UnwrapANGLE_line = UnwrapANGLE(Iy,:);
% unwrapping central horizontal line:

UnwrapANGLE_line(Ix:end) = unwrap(ANGLE(Iy,Ix:end));

UnwrapANGLE_line(1:(Ix-1)) = fliplr(unwrap(fliplr(ANGLE(Iy,1:(Ix-1)))));


ANGLE = unwrap(ANGLE ,2);


% 
% figure;
% plot(ANGLE(Iy,:))
% hold on
% plot(UnwrapANGLE(Iy,:),'r')


% adjusting other direction:
UnwrapANGLE = ANGLE - ones(size(ANGLE,1),1)*(ANGLE(Iy,:)-UnwrapANGLE_line);

UnwrapANGLE = UnwrapANGLE+PhaseRef-UnwrapANGLE(Iy,Ix);

% figure;
% subplot(1,2,1)
% surf(ANGLE)
% shading interp
% subplot(1,2,2)
% surf(UnwrapANGLE)
% shading interp




end

