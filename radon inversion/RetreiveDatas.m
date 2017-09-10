function datas = RetreiveDatas( raw , Ntrig, Nlines,MedElmtList)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% test to see o raw as the rigth number of coulumn :
if size(raw,2)~=Ntrig*Nlines
    datas = [];
    fprintf('wrong num of columns for input matrix')
    
else
    
    
fprintf('For current processed data : \r Nlines = %d , Ntrig = %d \r\n',Nlines, Ntrig)

[Isorted,Iposition] = sort(MedElmtList);

datas = reshape(raw,size(raw,1),Nlines, Ntrig);
datas(:,Isorted,:) = datas(:,Iposition,:);
datas = mean(datas,3);
end


end

