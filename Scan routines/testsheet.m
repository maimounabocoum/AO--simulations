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

imagesc(Datas)

figure; plot(Datas(:))












