%	Date: 2015
%	Author: Maïmouna Bocoum
%
%*********************************************************
%  function description:

% exemples:
% [D_out x_out y_out z_out] = ReduceDataSize(D_in ,'x',x_in,[xmin xmax],'y',y_in,[ymin ymax]...)
% [D_out x_out y_out] = ReduceDataSize(D_in x_in y_in,'xlim',[xmin xmax],)
% [D_out x_out] = ReduceDataSize(D_in x_in,'x',x_in,[xmin xmax]) : D_in is
% the data matrix to reduce, x_in to be reduced between the values xmin and
% xmax

function varargout = ReduceDataSize(varargin)


narginchk(4,10);

% number of input arguments sould be correct:
if mod(nargin -1,3)>0
   error(message('You are missing some argument. Make shure you have entered (Data,string,vector,[min max])')) 
end


DataMatrix = varargin{1};


Istr = intersect(1:(nargin -1),[1 4 7]);

% loop on axis vectors to extract:
for i = 1:length(Istr)

    Input_string = varargin{Istr(i)+1};
    
    if ischar(Input_string) == 0
    error(message('Please enter strings values to extract x,y and/or z')) 
    end
    
    switch  Input_string
        
        case 'x'

            v = varargin{Istr(i)+2};
            Xlim = varargin{Istr(i)+3};
            DataMatrix = DataMatrix(:,v>=Xlim(1)&v<=Xlim(2),:);
            v = v(v>=Xlim(1)&v<=Xlim(2));
                     
        case 'y'
            
            v = varargin{Istr(i)+2};
            Ylim = varargin{Istr(i)+3};
            DataMatrix = DataMatrix(v>=Ylim(1)&v<=Ylim(2),:,:);
            v = v(v>=Ylim(1)&v<=Ylim(2));
        case 'z'  
            
            v = varargin{Istr(i)+2};
            Zlim = varargin{Istr(i)+3};
            DataMatrix = DataMatrix(:,:,v>=Zlim(1)&v<=Zlim(2));
            v = v(v>=Zlim(1)&v<=Zlim(2));
    end
    
        varargout{i+1} = v;
    
end

%Nmatrix = size(varargin{1})

varargout{1} = DataMatrix;


end