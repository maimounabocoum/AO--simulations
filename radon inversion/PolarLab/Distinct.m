function [Number]=Distinct(Array)

%=======================================================================
% This function counts the number of distict and different ROWS in Array. It is 
% good in order to find how many repetitions we have in a created grid 
% (see example)
%
% Synopsis:  [Number]=Distinct(Array)
%
% Input:
%       Array - a 2-D array with possible repetitions of rows in it. 
% Output:
%       Number - the count of the distinct values
%
% Example:
%       [XC,YC]=Create_Grid('P',[20,20,sqrt(2)*pi],'');
%       disp(Distinct([XC(:),YC(:)]));
%
% Written by Michael Elad on Mrch 20th, 2005. 
%=======================================================================

M=size(Array,1);

Number=1;
for k=2:1:M,
    err=sum(abs(Array(1:k-1,:)-ones(k-1,1)*Array(k,:)),2);
    if min(err)>1e-8, 
        Number=Number+1;
    end;
end;

return;