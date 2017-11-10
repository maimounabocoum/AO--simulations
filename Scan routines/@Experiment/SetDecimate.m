function BoolActiveList = SetDecimate(obj,fx,BoolActiveList,type) 

Xm = obj.MyProbe.center(:,1);

switch type
    case 'sin'

BoolActiveList = repmat( (sin(2*pi*fx*Xm)>0) , 1, size(BoolActiveList,2) );
    case 'cos'
BoolActiveList = repmat( (cos(2*pi*fx*Xm)>0) , 1, size(BoolActiveList,2) );
end


end

