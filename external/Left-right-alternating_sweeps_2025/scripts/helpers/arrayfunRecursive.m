function [B, hasErrors] = arrayfunRecursive(A, fcn, targetClasses, defaultValue, catchErrors)
% Evaluate a function recursively on all applicable elements in array
if nargin < 3 || isempty(targetClasses), targetClasses = []; end
if nargin < 4 || isempty(defaultValue), defaultValue = []; end
if nargin < 5 || isempty(catchErrors), catchErrors = nargin >= 4; end

B = A;

% recursiveStruct = ~any(strcmp(targetClasses, "struct"));
% recursiveCell = ~any(strcmp(targetClasses, "cell"));

% gpuArray/strings are a special case, since they don't hold references to
% other data. So we treat them like basic data arrays.
AisObj = isobject(A) && ~isgpuarray(A) && ~isstring(A);

AmatchesTargetClass = any(strcmp(targetClasses, class(A)));

hasErrors = false;

% Is the input a struct/object array? If so, iterate through fields and
% process each as a separate array. *Unless* the input class matches the
% user-specified target class (in which case we evaluate the function on
% the data)
if (isstruct(A) || AisObj) && ~AmatchesTargetClass
    fds = string(fieldnames(A))';
    for fd = fds
        vals = {A.(fd)};
        % Call recursively on this sub-array. Wrapping it as a cell array
        % array, means we'll also get it back as a cell array
        [vals, hasErrorsTmp] = arrayfunRecursive(vals, fcn, targetClasses, defaultValue, catchErrors);
        [B.(fd)] = deal(vals{:});
        hasErrors = hasErrors | hasErrorsTmp;
    end
elseif iscell(A) && ~AmatchesTargetClass
    % Handle cell arrays by iterating through cells, calling recursively on
    % the contents of one cell at a time. As above, if "cell" is the
    % user-specified target class, we don't call recursively and instead
    % evaluate the user function on the array.
    for i = 1:numel(A)
        [B{i}, hasErrorsTmp] = arrayfunRecursive(A{i}, fcn, targetClasses, defaultValue, catchErrors);
        hasErrors = hasErrors | hasErrorsTmp;
    end
else
    % We reach here either (1) by default, if A is neither a struct/object 
    % nor a cell array, or (2) if class(A) matches a user-specified target
    % class.
    %
    % Either way, we have data of the target class. Evaluate the user
    % function. Optionally we can catch errors if they occur here, and
    % replace the output with a default value.
    if isempty(targetClasses) || AmatchesTargetClass
        if catchErrors
            try
                B = fcn(A);
            catch
                B = defaultValue;
                hasErrors = true;
            end
        else
            B = fcn(A);
        end
    end
end

end