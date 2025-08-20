function [paths, names] = filesearch_regexp(startDir, expression, recurse)
% filesearch_regexp searches for files in a directory that match a given regular expression.
%
% Usage:
%   [paths, names] = filesearch_regexp(startDir, expression, [norecurse])
%
% Inputs:
%   startDir   - Starting directory for the search
%   expression - Regular expression pattern to match file names
%   recurse  - Optional flag to control recursive search (default: 1 for recursive)
%
% Outputs:
%   paths - Cell array of paths to files matching the expression
%   names - Cell array of names of files matching the expression
%
% Example:
% searchPath='W:\HC_EEG'
% filepattern = '.*mid.*';
% [filepath,filename]=filesearch_regexp(searchPath, filepattern,0)
%
% 2024/02/06: fix return directory (e.g., mff folder from egi device)
% Author: Zhipeng Cao; zhipeng30@foxmail.com

if nargin < 3 || isempty(recurse)
    recurse = 1; % Default is to search recursively
end

% Check if the start directory exists
if ~exist(startDir, 'dir')
    error('Starting directory does not exist.');
end

% Internal function to perform the directory search
[paths, names] = dir_search(startDir, expression, recurse);

    function [paths, names] = dir_search(currDir, expression, recurse)
        paths = {};
        names = {};
        list_currDir = dir(currDir); % List all items in the current directory
        
        for u = 1:length(list_currDir)
            itemName = list_currDir(u).name;
            itemPath = fullfile(currDir, itemName); % Build full path for the item
            
            
            
            if recurse&& list_currDir(u).isdir && ~strcmp(itemName, '.') && ~strcmp(itemName, '..')
                
                % Recursively search in subdirectories if recurse is set
                [subPaths, subNames] = dir_search(itemPath, expression, recurse);
                paths = [paths, subPaths];
                names = [names, subNames];
                
            else
                % Check if the item matches the regular expression
                if ~isempty(regexpi(itemName, expression))
                    paths = [paths, {currDir}];
                    names = [names, {itemName}];
                end
            end
        end
    end
end
