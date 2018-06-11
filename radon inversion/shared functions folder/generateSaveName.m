function filenameOUT = generateSaveName( varargin )

% study the input vargin struture :
mod(nargin,2)
% check for an uneven number of input
if mod(nargin,2)~= 1 && nargin > 0
    
    filenameOUT = 'NoName' ;
    
else
    
 SubFolderName = varargin{1} ; 
 filenameOUT   = [SubFolderName,'\'] ;
 
 if nargin > 1
     % populate input list
     for Ninput = 1:(nargin - 1)/2

         InputName{Ninput}   = varargin{2*Ninput} ;
         InputValue{Ninput}  = num2str( varargin{2*Ninput+ 1} );

     end
     
     % put the field 'name in fist position ' :
     index = find(ismember( InputName , 'name')) ; 
     filenameOUT = strcat( filenameOUT, InputValue{index} );
     
     % remove element from cells
     InputValue(index) = [] ;
     InputName(index)  = []  ;
 
     for index = 1:length(InputName)
     filenameOUT = strcat( filenameOUT, '_' , InputName{index} ,'_', InputValue{index} );   
     end

     
 end
 
 
 

% get curent date


strdate = datestr( now,'_HH_MM_SS') ; 
strdate(4) = strrep(strdate(4),'_','h') ;  


filenameOUT = strcat( filenameOUT, strdate );


end


end

