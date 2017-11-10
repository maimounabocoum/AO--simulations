%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  main  program  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% maimouna bocoum 04-01-2017

Nmax = 5000 ;

Norm = repmat(1:Nmax,672,1);
Data = cumsum(S.Lines(:,1:Nmax),2)./Norm ;
imagesc(Data)
plot(Data(149,:))
hold on
plot(Data(500,:))

figure
plot(Data(:,1))
hold on
plot(smooth(Data(:,1),100))
hold on
plot(Data(:,Nmax))

%% 2D datas
% N = 2^11;
% Fmax = 1/(1000e-3);
% P = TF_t(N,Fmax);


figure; imagesc(Datas)
imagesc(imgaussfilt(Datas,3))


filter0 = exp(-(z-10).^2/(2)^2);
filter = repmat(filter0,1,size(Datas,2));
FDatas = fft(Datas(:)');
Ffilter = fft(filter);
Ffilter = ifftshift(imgaussfilt(fftshift(abs(Ffilter)),5));

figure; plot(fftshift(abs(FDatas))/max(abs(FDatas))) 
hold on 
plot(fftshift(abs(Ffilter)./max(abs(Ffilter))))

% FDatas(6) = 0 ;
% FDatas(28) = 0 ;
% FDatas(129) = 0 ;
% FDatas(203) = 0 ;
% 
% FDatas(end-6-1) = 0 ;
% FDatas(end-28-1) = 0 ;
% FDatas(end-129-1) = 0 ;
% FDatas(end-203-1) = 0 ;

FDatas = FDatas.*abs(Ffilter);



RDatas = ifft(FDatas,'symmetric');
RDatas = reshape(RDatas,[size(Datas,1),size(Datas,2)]);

figure; imagesc(RDatas)
%imagesc(imgaussfilt(RDatas,4))

%% data analysis for OP
    Hresconstruct = figure;
    set(Hresconstruct,'WindowStyle','docked');
    % plotting delay map
%     for i = 1:size(Delay,1)
%       Z_m(i,:) =   -Delay(i,:)*c*1e-6 ;
%     end

   [MconvX,MconvY] = meshgrid(-10:10,10:10);
    Mconv = exp(-MconvX.^2/(2*5)^2-MconvY.^2/(2*5)^2);
    %imagesc(imgaussfilt(Datas,3))
    A = imgaussfilt(Datas,5);
    B = conv2(Datas,Mconv,'same');
    MyImage = OP(A,Alphas,z,SampleRate*1e6,c) ;
    [I,z_out] = DataFiltering(MyImage) ;
%     Xm = (1:system.probe.NbElemts)*(0.2e-3) ;
    Xm = (1:192)*(0.2*1e-3) ;
    [theta,M0,X0,Z0] = EvalDelayLaw_shared(Xm,Delay,ActiveLIST,c);   
    
    Retroprojection_shared(I , Xm , z_out ,theta ,M0,Hresconstruct);
    % RetroProj_cleaned(Alphas,Datas,SampleR












