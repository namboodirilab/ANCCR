function fnames = findfiles(targetstr,directory,recurse,targetdir)
% Finds all files with the specified extension in the current directory and
% subdirectories (recursive). Returns a cell array with the fully specified
% file names.
%
%    files = findfiles( ext )
%       searches in the current directory and subdirectories.
%
%    files = findfiles( ext, directory)
%       starts search in directory specified.
%
%    files = findfiles( ext, directory, searchSubdirectories )
%       If searchSubdirectories == 0, the function will not descend
%       subdirectories. The default is that searchSubdirectories ~= 0.
%
% Note: if this function causes your recursion limit to be exceeded, then
% it is most probably trying to follow a symbolic link to a directory that
% has a symbolic link back to the directory it came from.
%
% Examples:
%
%    bmpfiles = findfiles('bmp')
%    dllfiles = findfiles('dll','C:\WinNT')
%    mfiles = findfiles('m',matlabroot)

% Copyright (C) 2001 Quester Tangent Corporation
% Author: Tony Christney, tchristney@questertangent.com
% $Id: findfiles.m 1.5 2001/05/25 01:54:17 tchristney Exp $

% Revision 1.6 2022/10/08 huijeong
% Now it searches targetstr in filename
%
% $Log: findfiles.m $
% Revision 1.5  2001/05/25 01:54:17  tchristney
% Fixed bug where numMatches was not updated after recursive call, causing
% potential overwrites in the fnames cell array.
%
% Revision 1.4  2001/05/24 20:46:09  tchristney
% Recursion on/off support courtesy Brett Shoelson.
% <Brett.Shoelson@joslin.harvard.edu>
%
% Revision 1.2  2001/03/16 21:28:37  tchristney
% added some comments for external release, i.e. examples.
%
% Revision 1.1  2001/03/07 22:38:52  tchristney
% Initial revision
%

if nargin == 1
    directory = cd;
    recurse = (1==1);
    targetdir = nan;
elseif nargin == 2
    oldDir = cd;
    cd(directory);
    directory = cd;
    cd(oldDir);
    recurse = (1==1);
    targetdir = nan;
elseif nargin >= 3
    oldDir = cd;
    cd(directory);
    directory = cd;
    cd(oldDir);
    if nargin ==3
        targetdir = nan;
    end
end

d = dir(directory);

fnames = {};
numMatches = 0;
try
for i=1:length(d)
    % look for occurences of targetstr in the file name
    if isnan(targetdir)
        extIndices = findstr(targetstr,d(i).name);
    else
        extIndices = findstr(targetstr,d(i).name) & findstr(targetdir,d(i).folder);
    end
    
    % if the file is not a directory, and the file has at least one occurence
    if ~d(i).isdir & ~isempty(extIndices)
        
        % then if the last occurrence is at the end of the file name,
        % add the file name to the list.
        %         if length(d(i).name) == (extIndices(length(extIndices)) + length(targetstr))
        numMatches = numMatches + 1;
        fnames{numMatches} = fullfile(directory,d(i).name);
        %         end
        % otherwise, descend directories appropriately.
        % note that this could result in a recursion limit error if it tries to
        % follow symbolic links that loop back on themselves...
    elseif recurse & d(i).isdir & ~strcmp(d(i).name,'.') & ~strcmp(d(i).name,'..')
        fnames = [fnames findfiles(targetstr,fullfile(directory,d(i).name),recurse,targetdir)];
        numMatches = length(fnames);
    end
end
catch
   i 
end
fnames = fnames(:);





