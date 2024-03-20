close all;
clear all;
clc

list_factory = fieldnames(get(groot, 'factory'));
index_interpreter = find(contains(list_factory, 'Interpreter'));

for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)}, 'factory', 'default');
    set(groot, default_name, 'latex');
end

filepath = 'BeamShapeEvolution_2TSOA.mat';
timeduringpulse = load(filepath, "timeduringpulse");
timeduringpulse = struct2cell(timeduringpulse);
timeduringpulse = cell2mat(timeduringpulse);

%% plot slow and fast axis - and respective fit
wx = load(filepath, "wx");
wx = struct2cell(wx);
wx = cell2mat(wx);

wy = load(filepath, "wy");
wy = struct2cell(wy);
wy = cell2mat(wy);

f=polyfit(timeduringpulse,wx,4) ;
wx_fit = polyval(f,timeduringpulse); 
g=polyfit(timeduringpulse,wy,12) ;
wy_fit = polyval(g,timeduringpulse); 


figure(1) ; clf ; 
plot(timeduringpulse,wx,'o','color','red') ; hold on 
plot(timeduringpulse,wx_fit,'color','red') ; hold on
plot(timeduringpulse,wy_fit,'color','blue') ; hold on
plot(timeduringpulse,wy,'o','color','blue') ;

xlabel('time (us)')
ylabel('waist(mm)')
grid on
%% Slow axis y

wx = wx_fit ;
wy = wy_fit ;
 
lambda = 828e-9;

f5 = 2.75e-3;
f6 = 12.7e-3;
D = 85e-2;
delta = (0:0.1:20).*1e-3;
w0 = 1.5e-6;
s = 640e-6 + f5 + delta;
Zr = pi*w0^2/lambda;

s_ = 1./(1./(s+Zr^2./(s+f5))+1/f5); d0 = s_;
w0_ = w0./(sqrt((1-s./f5).^2+(Zr./f5).^2));
Zr_ = pi*w0_.^2./lambda;

s__ = 1./(1./(f6+Zr_.^2./(f6+f6))+1./f6);
w0__ = w0_./(sqrt((1-f6./f6)^2+(Zr_./f6).^2));

z= D - s__;
w2 = w0__.*sqrt(1+(lambda*z./pi./w0__.^2).^2);

plot(delta*1e3,w2*1e3)

p = polyfit(delta, w2,5);

for i = 1:length(wy)
% Valeur de x pour laquelle on veut estimer y
x_value = wy(i)*1e-3;

% Évaluation du polynôme pour x_value
delta_values(i) = polyval(p, x_value);

end

plot(delta_values)

z_ref = delta_values + 640e-6;
R_ref = z_ref.*(1+(pi*w0^2/lambda./z_ref).^2);
w_ref = w0.*sqrt(1+(z.*lambda./pi./w0^2).^2);
plot(R_ref)
plot(w_ref)


% Définition de la grille spatiale
Nx = 512; % Nombre de points en x
Ny = 512; % Nombre de points en y
Lx = 100e-2; % Taille de la grille en x en m
Ly = 100e-2; % Taille de la grille en y en m
dx = Lx / Nx; % Pas d'échantillonnage en x
dy = Ly / Ny; % Pas d'échantillonnage en y
x = linspace(-Lx/2, Lx/2, Nx); % Grille spatiale en x
y = linspace(-Ly/2, Ly/2, Ny); % Grille spatiale en y
[X, Y] = meshgrid(x, y);


%% Retreival of the input Field

for i = 1:length(wy)

    Profile{i}  = exp(-( (X.^2)/(wx(i)^2) + Y.^2/(wy(i)^2) ))   ;
    phaseLaw{i} = exp(1i*(X.^2+Y.^2)./2/R_ref(i))               ;

    

figure(1); clf ;
subplot(121)
imagesc(x*1e3,y*1e3,Profile{i})
xlabel('x(mm)')
ylabel('y(mm)')
subplot(122)
imagesc(x*1e3,y*1e3,angle(phaseLaw{i}))
xlabel('x(mm)')
ylabel('y(mm)')
title(num2str(i))
end

%% field cross-correlation
for i = 1:length(wy)
    Eref = Profile{100}.*phaseLaw{100};
    E  = Profile{i}.*phaseLaw{i} ;
    
    g1(i) = ( sum( E(:).*conj(Eref(:))) )/sqrt((sum( E(:).*conj(E(:))) )*(sum( Eref(:).*conj(Eref(:))) )) ;

end


figure(2)
plot(abs(g1))
%%














