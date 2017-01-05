function [positions, amp] = phantom(N)

x_size = 20/1000;
y_size = 20/1000;
z_size = 20/1000;
z_start = 30/1000;

x = (rand(N,1)-0.5)*x_size;
y = (rand(N,1)-0.5)*y_size;
z = rand(N,1)*z_size + z_start;

% x = (ones(N,1)-0.5)*x_size;
% y = (ones(N,1)-0.5)*y_size;
% z = ones(N,1)*z_size + z_start;

amp = randn(N,1);

% r = 1.5/1000;
% xc = 0/1000;
% zc = 10/1000 + z_start;
% 
% inside = ( ((x-xc).^2 + (z-zc).^2) < r^2);
% amp = amp.* (1-inside);
% 
% dz= z_size/10;
% for i=N-9:N
%     x(i) = 0;%-15/1000;
%     y(i) = 0;
%     z(i) = z_start + (i-N+9)*dz;
%     amp(i) = 100;
% end

positions = [x y z];

end



