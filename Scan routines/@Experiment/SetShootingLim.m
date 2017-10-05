function obj = SetShootingLim(obj)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% reindexation of X0, X1 for 0 to Lmax of the probe
X0 = obj.param.X0 + (1/2)*obj.param.N_elements*obj.param.width;
X1 = obj.param.X1 + (1/2)*obj.param.N_elements*obj.param.width;
      
ElmtBorns   = [min(obj.param.N_elements,max(1,round(X0/obj.param.width))),...
               max(1,min(obj.param.N_elements,round(X1/obj.param.width)))];
ElmtBorns   = sort(ElmtBorns); % in case X0 and X1 are mixed up

I = 1:obj.param.N_elements ;
I = setdiff(I, ElmtBorns(1):ElmtBorns(2)) ; % remove element not in shooting defined window


% set excluded actuators to false :
obj.BoolActiveList(I,:) = false ;

end

