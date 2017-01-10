%==========================================================================
%                           This script creates the figures for the Polar_Paper.TeX
%==========================================================================

% Figure 1 - The polar coordinates system in the square [-pi,pi]

figure(1); clf; 
N=8;
plot(0,0,'.w'); hold on;
for k=0:1:N,
    t=0:pi/1000:2*pi;
    r=k*pi/N;
    plot(r*cos(t),r*sin(t),'b:');
end;
Points=0:pi/(2*N):2*pi-eps;
for k=1:1:length(Points), 
    plot([0,pi*cos(Points(k))],[0,pi*sin(Points(k))],'b:');
end;
hold on;
axis equal; axis([-4 4 -4 4]);
[GridX,GridY]=Create_Grid('P',[9,32,9*pi/8],'');
h=plot(GridX,GridY,'ob');
set(h,'MarkerFaceColor','b');
xlabel('\xi_x'); ylabel('\xi_y'); 
plot([-pi,pi,pi,-pi,-pi],[-pi,-pi,pi,pi,-pi],'b'); 
print -depsc2 Figure001.eps;

% Figure 2 - The Cartesian-to-Polar and Polar-to-Cartesian in USFFT

figure(1); clf; 
subplot(1,2,1); 
plot([-pi,pi,pi,-pi,-pi],[-pi,-pi,pi,pi,-pi],'g'); 
hold on;
xlabel('\xi_x'); ylabel('\xi_y'); title('Cartesian-to-Polar');
Create_Grid('C',[32,32,[-pi,pi],[-pi,pi]],'.g'); axis([-4 4 -4 4]);
[x,y]=Create_Grid('P',[9,16,9*pi/8],''); 
h=plot(x,y,'ob'); 
set(h,'MarkerFaceColor','b')

subplot(1,2,2); 
[x,y]=Create_Grid('P',[17,32,17*pi/16],''); 
h=plot(x,y,'.g'); 
hold on; 
plot([-pi,pi,pi,-pi,-pi],[-pi,-pi,pi,pi,-pi],'g'); 
xlabel('\xi_x'); ylabel('\xi_y'); title('Polar-to-Cartesian');
[x,y]=Create_Grid('C',[8,8,[-pi,pi],[-pi,pi]],'');
h=plot(x,y,'ob'); 
set(h,'MarkerFaceColor','b')
axis image; axis([-4 4 -4 4]);
print -depsc2 Figure002.eps;

% Figure 3 - The pseudo-polar grid

figure(1); clf; 

Create_Grid('R',[9,16,9*pi/8],'.b'); 
hold on;
[x,y]=Create_Grid('R',[9,16,9*pi/8],'');
h=plot(x,y,'ob'); 
set(h,'MarkerFaceColor','b')

xlabel('\xi_x'); ylabel('\xi_y');
plot(1*[-pi,pi,pi,-pi,-pi]/8,1*[-pi,-pi,pi,pi,-pi]/8,'b'); 
plot(2*[-pi,pi,pi,-pi,-pi]/8,2*[-pi,-pi,pi,pi,-pi]/8,'b'); 
plot(3*[-pi,pi,pi,-pi,-pi]/8,3*[-pi,-pi,pi,pi,-pi]/8,'b'); 
plot(4*[-pi,pi,pi,-pi,-pi]/8,4*[-pi,-pi,pi,pi,-pi]/8,'b'); 
plot(5*[-pi,pi,pi,-pi,-pi]/8,5*[-pi,-pi,pi,pi,-pi]/8,'b'); 
plot(6*[-pi,pi,pi,-pi,-pi]/8,6*[-pi,-pi,pi,pi,-pi]/8,'b'); 
plot(7*[-pi,pi,pi,-pi,-pi]/8,7*[-pi,-pi,pi,pi,-pi]/8,'b'); 
plot(8*[-pi,pi,pi,-pi,-pi]/8,8*[-pi,-pi,pi,pi,-pi]/8,'b'); 
axis([-4 4 -4 4]);
plot([-pi,pi],[-pi,pi],'b');
plot(2*[-pi,pi]/4,[-pi,pi],'b');
plot(0*[-pi,pi]/4,[-pi,pi],'b');
plot(-2*[-pi,pi]/4,[-pi,pi],'b');
plot([pi,-pi],[-pi,pi],'b');
plot([pi,-pi],2*[-pi,pi]/4,'b');
plot([pi,-pi],0*[-pi,pi]/4,'b');
plot([pi,-pi],-2*[-pi,pi]/4,'b');
print -depsc2 Figure003.eps;

% Figure 5 - The Pseudo-Polar Grid - separation to BV and BH

figure(1); clf; 
plot(0,0,'.w'); hold on;
xlabel('\xi_x'); ylabel('\xi_y');
for k=1:1:8,
    plot(k*[-pi,pi,pi,-pi,-pi]/8,k*[-pi,-pi,pi,pi,-pi]/8,'b'); 
end;
for k=-1+1/4:1/4:1,
    plot([-pi pi],k*[-pi pi],'b');
    plot(-k*[-pi pi],[-pi pi],'b');
end;

N=8;
for ll=-N:1:N-1, % the BV
    y=pi*ll/N;
    for m=-N/2:1:N/2-1,
        x=2*pi*ll*m/N^2;
        h=plot(x,y,'ob');
        set(h,'MarkerFaceColor','b')
    end;
end;
for ll=-N:1:N-1, % the BH
    x=pi*ll/N;
    for m=-N/2+1:1:N/2,
        y=2*pi*ll*m/N^2;
        h=plot(x,y,'ob');
        set(h,'MarkerFaceColor','w')
    end;
end;
axis equal; axis([-4 4 -4 4]);
print -depsc2 Figure005.eps

% Figure 9 - The Polar Grid - separation to BV and BH

figure(2); clf; 
N=8;
plot(0,0,'.w'); hold on;
xlabel('\xi_x'); ylabel('\xi_y');
for k=0:1:N,
    t=0:pi/1000:2*pi;
    r=k*pi/N;
    plot(r*cos(t),r*sin(t),'b');
end;
Points=0:pi/(2*N):2*pi-eps;
for k=1:1:length(Points), 
    plot([0,pi*cos(Points(k))],[0,pi*sin(Points(k))],'b');
end;
axis equal; axis([-4 4 -4 4]);

for r=-pi:pi/N:pi-pi/N, %BV
    for t=-pi/4:pi/(2*N):pi/4-pi/(2*N),
        y=r*cos(t); 
        x=r*sin(t);
        h=plot(x,y,'ob');
        set(h,'MarkerFaceColor','b')
    end;
end;
for r=-pi:pi/N:pi-pi/N, %BH
    for t=pi/4:pi/(2*N):3*pi/4-pi/(2*N),
        y=r*cos(t); 
        x=r*sin(t);
        h=plot(x,y,'ob');
        set(h,'MarkerFaceColor','w')
    end;
end;
print -depsc2 Figure009.eps

% Figure 10 - Treatment of the rays in the PP to Polar
figure(1); clf; 
plot(0,0,'.w'); hold on;
xlabel('\xi_x'); ylabel('\xi_y');
for k=1:1:8,
    plot(k*[-pi,pi,pi,-pi,-pi]/8,k*[-pi,-pi,pi,pi,-pi]/8,'b:'); 
end;
for k=-1+1/4:1/4:1,
    plot([-pi pi],k*[-pi pi],'b:');
    plot(-k*[-pi pi],[-pi pi],'b:');
end;
axis image; axis([-4 4 -4 4]);

for theta=-pi/4:pi/16:pi/4-pi/16,
     h=plot([-pi,pi],[-tan(theta),tan(theta)]*pi,'k');   
     set(h,'LineWidth',1.5);
end;
for theta=-pi/4:pi/16:pi/4-pi/16,
     h=plot([tan(theta),-tan(theta)]*pi,[-pi,pi],'k');   
     set(h,'LineWidth',1.5);
end;
print -depsc2 Figure010.eps

% Figure 11 - the first interpolation stage
figure(1); clf; 
plot(0,0,'.w'); hold on;
xlabel('\xi_x'); ylabel('\xi_y');
for k=1:1:8,
    plot(k*[-pi,pi,pi,-pi,-pi]/8,k*[-pi,-pi,pi,pi,-pi]/8,'b:'); 
end;
for k=-1+1/4:1/4:1,
    plot([-pi pi],k*[-pi pi],'b:');
    plot(-k*[-pi pi],[-pi pi],'b:');
end;
axis image; axis([-4 4 -4 4]);

for theta=-pi/4:pi/16:pi/4-pi/16,
     h=plot([-pi,pi],[-tan(theta),tan(theta)]*pi,'b');   
     set(h,'LineWidth',1.5);
end;
for theta=-pi/4:pi/16:pi/4-pi/16,
     h=plot([tan(theta),-tan(theta)]*pi,[-pi,pi],'b');   
     set(h,'LineWidth',1.5);
end;

N=8;
ll=N-2; % one row from the BV
for m=-N/2:1:N/2-1,
     x=2*pi*ll*m/N^2;
     h=plot(x,y,'.b');
     set(h,'MarkerSize',35);
end;
y=pi*ll/N;
for theta=-pi/4:pi/16:pi/4-pi/16,
    x=y*tan(theta);
     h=plot(x,y,'sb');
     set(h,'MarkerFaceColor','w');
end;
print -depsc2 Figure011.eps

% Figure 12 - Treatment of circles 
figure(1); clf; 
plot(0,0,'.w'); hold on;
xlabel('\xi_x'); ylabel('\xi_y');
for k=1:1:8,
    plot(k*[-pi,pi,pi,-pi,-pi]/8,k*[-pi,-pi,pi,pi,-pi]/8,'b:'); 
end;
axis image; axis([-4 4 -4 4]);

for theta=-pi/4:pi/16:pi/4-pi/16,
     h=plot([-pi,pi],[-tan(theta),tan(theta)]*pi,'b:');   
end;
for theta=-pi/4:pi/16:pi/4-pi/16,
     h=plot([tan(theta),-tan(theta)]*pi,[-pi,pi],'b:');   
end;
N=8;
for k=0:1:N,
    t=0:pi/1000:2*pi;
    r=k*pi/N;
    plot(r*cos(t),r*sin(t),'b');
end;
print -depsc2 Figure012.eps

% Figure 13 - Treatment of circles - the second interpolation 
figure(1); clf; 
plot(0,0,'.w'); hold on;
xlabel('\xi_x'); ylabel('\xi_y');
for k=1:1:8,
    plot(k*[-pi,pi,pi,-pi,-pi]/8,k*[-pi,-pi,pi,pi,-pi]/8,'b:'); 
end;
axis image; axis([-4 4 -4 4]);
for theta=-pi/4:pi/16:pi/4-pi/16,
     h=plot([-pi,pi],[-tan(theta),tan(theta)]*pi,'b:');   
end;
for theta=-pi/4:pi/16:pi/4-pi/16,
     h=plot([tan(theta),-tan(theta)]*pi,[-pi,pi],'b:');   
end;
N=8;
for k=0:1:N,
    t=0:pi/1000:2*pi;
    r=k*pi/N;
    plot(r*cos(t),r*sin(t),'b');
end;
for ll=-N:1:N-1, % the BV
    y=pi*ll/N;
    x=y*tan(-pi/4+pi/16);
    h=plot(x,y,'.b');
    set(h,'MarkerSize',35);
end;
for ll=-N:1:N, % the BV
    y=pi*ll/N;
    x=y*tan(-pi/4+pi/16);
    x=x*cos(pi/4-pi/16);
    y=y*cos(pi/4-pi/16);
    h=plot(x,y,'bs');
    set(h,'MarkerFaceColor','w');    
end;
print -depsc2 Figure013.eps

% Figure 16 - Disk-Band-Limited signals 
figure(1); clf; 
plot(0,0,'.w'); hold on;
xlabel('\xi_x'); ylabel('\xi_y');
h=plot([-pi -pi pi pi -pi],[-pi pi pi -pi -pi],'r');
set(h,'LineWidth',2);
axis image;
axis([-4 4 -4 4]);
h=plot([-pi -pi pi pi -pi]/sqrt(2),[-pi pi pi -pi -pi]/sqrt(2),'r');
set(h,'LineWidth',2);
t=0:pi/100:2*pi,
h=plot(pi*sin(t),pi*cos(t),'b');
set(h,'LineWidth',2);
grid on;
print -depsc2 Figure016.eps

% Figure 21 - approximation error as a function of the oversampling factor
Comparison_1_USFT_vs_Polar(64,5);
load Comparison_1_Results.mat
figure(1); clf;
subplot(2,1,1); 
semilogy((1:1:20).^2,Err2(:,1),'--'); hold on;
semilogy((1:1:20)*1,Err3(:,1));
semilogy((1:1:20)*2,Err4(:,1));
semilogy((1:1:20)*3,Err5(:,1));
semilogy((1:1:20)*4,Err6(:,1));
semilogy((1:1:20)*5,Err7(:,1));
semilogy((1:1:20)*6,Err8(:,1));
axis([0 80 1e-7 1]);
gtext('USFFT');
gtext('Fast Polar with P=1,2, ... , 6');
title('L_1 norm approx. error');
xlabel('Over-Sampling factor');
subplot(2,1,2);
semilogy((1:1:20).^2,Err2(:,2),'--'); hold on;
semilogy((1:1:20)*1,Err3(:,2));
semilogy((1:1:20)*2,Err4(:,2));
semilogy((1:1:20)*3,Err5(:,2));
semilogy((1:1:20)*4,Err6(:,2));
semilogy((1:1:20)*5,Err7(:,2));
semilogy((1:1:20)*6,Err8(:,2));
axis([0 80 1e-7 1]);
gtext('USFFT');
gtext('Fast Polar with P=1,2, ... , 6');
title('L_2 norm (square-root MSE) approx. error');
xlabel('Over-Sampling factor');
figure(1); print -depsc2 Figure021.eps

% Figure 22
Comparison_2a_USFT_vs_Polar(64,4,20,5,0.01);
figure(1); print -depsc2 Figure022a.eps
figure(2); print -depsc2 Figure022b.eps

% Figure 23
Comparison_3_USFT_vs_Polar(16,4,20,3,0.02);
figure(2); print -depsc2 Figure023.eps

% Figures 24-29
Comparison_2_USFT_vs_Polar(16,4,20);
figure(1); print -depsc2 Figure024.eps
figure(2); print -depsc2 Figure025.eps
figure(3); print -depsc2 Figure026.eps
figure(4); print -depsc2 Figure027.eps
figure(5); print -depsc2 Figure028.eps
figure(6); print -depsc2 Figure029.eps

% Figure 32 - the bound on the error
C=2/((2*pi)^4*3);
for P=1:1:5,
    S=1:1:80/P,
    E=(2*pi)^4*C*(2./S.^4+1/P^4);
    if P==1, 
        figure(1); clf; semilogy(P*S,E,'r'); hold on;
    else,
        semilogy(P*S,E,'r'); 
    end;
end;
grid on;
xlabel('Over-sampling factor [S\cdot P]');
ylabel('L_1 Average Error Bound');
gtext('P=1');
gtext('P=2');
gtext('P=3');
gtext('P=4');
gtext('P=5');
print -depsc2 Figure032.eps

% Figure 33 - the effect of the corners
Comparison_4_Polar_Corners;
print -depsc2 Figure033.eps















