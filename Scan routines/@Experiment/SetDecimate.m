function BoolActiveList = SetDecimate(obj,decimation,BoolActiveList,type) 

Xm = obj.MyProbe.center(:,1);

switch type
    case 'sin'
BoolActiveList = repmat( (sin(2*pi*decimation*Xm)>0) , 1, size(BoolActiveList,2) );
    case 'cos'
BoolActiveList = repmat( (cos(2*pi*decimation*Xm)>0) , 1, size(BoolActiveList,2) );
end

%% old decimate version (to be erased)
if 1 == 0
%function to set the decimate of the probe where BoolActiveList is matrix 
% this function act on all line of that matrix
N_elements = size(BoolActiveList,1) ;

% calculate period from decimate:
period_pixels = min(N_elements + 1 ,N_elements/decimation ) ;

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
period_pixels = round(  period_pixels  ) ;

% center   = mean(ElmtBorns) + floor(decimation/2) ; % in case X0 and X1 are mixed up

Imod = mod( (1:N_elements) - ElmtBorns(1) , 2*period_pixels ) ;

                       
BoolActiveList( Imod < period_pixels , : )  = true;
BoolActiveList( Imod >= period_pixels , : ) = false ;
                    
end

end

