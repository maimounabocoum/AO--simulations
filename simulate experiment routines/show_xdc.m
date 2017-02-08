function res = show_xdc (Th)
% display for a linear transducer 
% Do it for the rectangular elements
colormap(cool(128));
data = xdc_get(Th,'rect');
%data = xdc_get(Th,'focus');
% M : number of rectangles
[N,M]=size(data);
% Do the actual display
for i=1:M
x=[data(11,i), data(20,i); data(14,i), data(17,i)]*1000;
y=[data(12,i), data(21,i); data(15,i), data(18,i)]*1000;
z=[data(13,i), data(22,i); data(16,i), data(19,i)]*1000;
c=data(5,i)*ones(2,2);
hold on
surf(x,y,z,c)
%Z(i) = data(13,i) ;
end
%hold on
%plot(Z)
% Put som axis legends on
Hc = colorbar;
view(3)
xlabel('x [mm]')
ylabel('y [mm]')
zlabel('z [mm]')
grid
axis('image')
hold off

end

