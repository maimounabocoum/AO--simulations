
 x_phantom = CurrentExperiement.MySimulationBox.x ;
 y_phantom = CurrentExperiement.MySimulationBox.y ;
 z_phantom = CurrentExperiement.MySimulationBox.z ;
 %[MyTansmission,R,zR] = CurrentExperiement.ShowPhantom(param.angles);
 
  GG = FourierMap(CurrentExperiement.ScanParam , param.nuX0 ,param.nuZ0 );
  GG = GG.DefineImage( MyTansmission );
  
  figure;plot(GG.In)
  
  %% TEST ARBITRARY NUMBER OF FFT
  
  s = rand(1,151);
  Fs = 1e3;
  F = TF_t(151,Fs);
  stf = F.fourier(s);
 
  figure;
  subplot(121)
  plot(F.t,s)
  subplot(122)
  plot(F.f,abs(stf))
  
  
  
  
  %%