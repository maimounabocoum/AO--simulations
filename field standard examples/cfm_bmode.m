%  Make a combined CFM and B-mode image from
%  simulated data
%
%  Version 1.0, 3/4-97, JAJ

%  Do the gray scale data first

make_image

%  Then do the cfm image

cfm_image

%  Combine the two images
%  Ndist is the number of samples between two velocity estimates

Ndist_new=D*Ndist*fn_bmode/fs;  %  New conversion factor between sampling frequencies

[Nsamples,Nlines]=size(new_env)
[Nest,Nlest]=size(new_est);

for line=1:Nlines
line
  for sample=1:Nsamples
    if (sample/Ndist_new+1 < Nest) 
      if (new_est(sample/Ndist_new+1,line) > 0.25*64)
        comb_image(sample,line)=new_est(sample/Ndist_new+1,line) + 64;
      else
        comb_image(sample,line)=new_env(sample,line)/2;
        end
    else
      comb_image(sample,line)=new_env(sample,line)/2;
      end
    end
  end
  
%  Make a color map for the combined image

map=[[0:63 0:63]; [0:63 zeros(1,64)]; [0:63 zeros(1,64)]]/63;
colormap(map')

%  Display the combined result

clg
dx=40/500;
image( ((1:Nlines)-Nlines/2)*dx,((1:Nsamples)/fn_bmode)*c/2*1000,comb_image)
xlabel('Lateral distance [mm]')
ylabel('Axial distance [mm]')
colormap(map')
%brighten(-0.35)
axis([-20 20 30 90])
axis('image')

print -depsc cfm_image.eps



