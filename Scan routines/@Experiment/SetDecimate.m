function BoolActiveList = SetDecimate(obj,decimation,BoolActiveList) 

%function to set the decimate of the probe where BoolActiveList is matrix 
% this function act on all line of that matrix

N_elements = size(BoolActiveList,1) ;

% periode : decimate should be a entire value
% center : index of the center for the that decimation

% reindexation of X0, X1 for 0 to Lmax of the probe
X0 = obj.param.X0 + (1/2)*obj.param.N_elements*obj.param.width;
X1 = obj.param.X1 + (1/2)*obj.param.N_elements*obj.param.width;



ElmtBorns   = [min(obj.param.N_elements,max(1,round(X0/obj.param.width))),...
               max(1,min(obj.param.N_elements,round(X1/obj.param.width)))];

% current window length in units
L = length(ElmtBorns(1):ElmtBorns(2)) ;
% decimation in units :
%decimation = round(  L/decimation  ) ;
decimation = round(  decimation  ) ;

% center   = mean(ElmtBorns) + floor(decimation/2) ; % in case X0 and X1 are mixed up

Imod = mod( (1:N_elements) - ElmtBorns(1) , 2*decimation ) ;

                       
BoolActiveList( Imod < decimation , : )  = true;
BoolActiveList( Imod >= decimation , : ) = false ;
                    


end

