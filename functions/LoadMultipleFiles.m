function [Filename,Foldername,Nfiles] = LoadMultipleFiles( folder )


% opens default folder :
if nargin == 0
[Filename,Foldername] = uigetfile('','MultiSelect','on');
Nfiles = length(Filename);
if Nfiles==1
   Filename = {Filename};
end

else
    
[Filename,Foldername] = uigetfile(folder,'MultiSelect','on');
Nfiles = length(Filename);
if Nfiles==1
   Filename = {Filename};
end    
    
end


end

