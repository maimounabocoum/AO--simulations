function folderOUT = generateSubFolderName( inputfolder )

% check if current folder exist :
%    0 — name does not exist.
%    1 — name is a variable in the workspace.
%    2 — name is a file with extension .m, .mlx, or .mlapp, or name is the name of a file with a non-registered file extension (.mat, .fig, .txt).
%    3 — name is a MEX-file on your MATLAB® search path.
%    4 — name is a Simulink® model or library file on your MATLAB search path.
%    5 — name is a built-in MATLAB function. This does not include classes.
%    6 — name is a P-code file on your MATLAB search path.
%    7 — name is a folder.
%    8 — name is a class. (exist returns 0 for Java classes if you start MATLAB with the -nojvm option.)    0 — name does not exist.
    
IsFolder = exist(inputfolder) ;

if IsFolder ~= 7
inputfolder =     uigetdir('') ;
end


% get curent date
Date = datestr(now,'yyyy-mm-dd') ;

% concatenate strings : 
folderOUT = strcat(inputfolder,'\',Date);

if ~exist(folderOUT)
IsCreated = mkdir(folderOUT) ;
switch IsCreated
    case 'true'
        fprinf('folder created successfully') ;
    case 'false'
        fprinf('folder not created') ;
end



end

