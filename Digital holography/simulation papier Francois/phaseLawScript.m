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

wx = load(filepath, "wx");
wx = struct2cell(wx);
wx = cell2mat(wx);

wy = load(filepath, "wy");
wy = struct2cell(wy);
wy = cell2mat(wy);

%% Slow axis y
 
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
Lx = 20e-2; % Taille de la grille en x en m
Ly = 20e-2; % Taille de la grille en y en m
dx = Lx / Nx; % Pas d'échantillonnage en x
dy = Ly / Ny; % Pas d'échantillonnage en y
x = linspace(-Lx/2, Lx/2, Nx); % Grille spatiale en x
y = linspace(-Ly/2, Ly/2, Ny); % Grille spatiale en y
[X, Y] = meshgrid(x, y);

for i = 1:length(wy)
    phaseLaw{i} = (exp(1i*(X.^2+Y.^2)./2/R_ref(i)))-mean2(exp(1i*(X.^2+Y.^2)./2/R_ref(i)));
end

%n0 = 10;
% for i = 0:length(phaseLaw)-1 - n0
%     correlation_ij=[];
%     for j = n0:length(phaseLaw)
%         if (i+j<length(phaseLaw))
%         correlation_ij = [correlation_ij...
%             abs(sum(sum(conj(phaseLaw{i+j}).*phaseLaw{j}))).^2/(sum(sum(abs(phaseLaw{i+j}).^2))*sum(sum(abs(phaseLaw{j}).^2)))];%corr2(phaseLaw{j},phaseLaw{i+j})
%         end
%     end 
%     xCorr(i+1) = mean(correlation_ij)
% end

for i = 1:length(phaseLaw)
    xCorr(i) = corr2(angle(phaseLaw{i}),angle(phaseLaw{100}));
    xCorr_g1(i) = abs(sum(sum(conj(phaseLaw{100}).*phaseLaw{i}))).^2/(sum(sum(abs(phaseLaw{100}).^2))*sum(sum(abs(phaseLaw{i}).^2)));
end



figure(1)
clf
plot(xCorr,'-*')
hold on
plot(xCorr_g1,'-*')
hold on
plot(sqrt(xCorr_g1),'-*')
hold on
% plot(xCorr_g2,'-*')
% hold on
plot(31:length(xCorr_exp)+30,xCorr_exp,'-*')
xlabel("Time [$\mu$s]")
grid on
grid minor
legend("Corrélation des lois de phase incidentes","$|g^{(1)}|^2$","$g^{(1)}$","Corrélation des speckles")
set(gca, 'FontSize', 14);