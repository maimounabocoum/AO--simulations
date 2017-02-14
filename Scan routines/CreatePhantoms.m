%%%phantom generation _ this is a test sheet and is meant to be erased
%% maimouna bocoum 04-01-2017
clearvars ;

addpath('..\Field_II')
field_init(0);


parameters;

Ph = Phantom(param.phantom.Positions,param.phantom.Sizes,param.phantom.Types);
FBox = AO_FieldBox(param.Xrange,param.Yrange,param.Zrange,param.Nx,param.Ny,param.Nz);
I_abs = Ph.CalculatePhantom(FBox.x,FBox.y,FBox.z) ;


figure ;
T = reshape(I_abs,[param.Ny,param.Nx,param.Nz]);
imagesc(FBox.x*1e3,FBox.z*1e3,squeeze(T)')
ylabel('z')
xlabel('x')
colorbar
 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% End Program - Free memory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

field_end;