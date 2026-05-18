function dirs = findDirs(rootDir, ignoreList, useRecursion)
% Generates a recursive list of all dirs beneath the specified dir,
% excluding any folders named '.git'.
if nargin < 2 || isempty(ignoreList), ignoreList = {}; end
if nargin < 3 || isempty(useRecursion), useRecursion = false; end
dirs = findDirsLocal(rootDir, ignoreList, useRecursion);

% Return in same 
if isstring(rootDir)
    dirs = string(dirs);
else
    dirs = cellstr(dirs);
end

end

function dirs = findDirsLocal(parentDir, ignoreList, useRecursion)
dirs = {};
% Get directory contents
dirContents = dir(parentDir);
% Loop through the contents
for k = 1:length(dirContents)
    % Ignore '.' and '..' and any hidden '.git' folders
    if dirContents(k).isdir && ~strcmp(dirContents(k).name, '.') && ~strcmp(dirContents(k).name, '..') ...
            && ~any(strcmp(dirContents(k).name, ignoreList))
        % Add the current directory to the list
        folderPath = fullfile(parentDir, dirContents(k).name);
        dirs = [dirs; {folderPath}];
        if useRecursion
            % Recursively get folders in the current directory
            dirs = [dirs; findDirsLocal(folderPath, ignoreList, useRecursion)];
        end
    end
end
end