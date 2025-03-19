function batchCompute()
% batchCompute  Processes all .fig files in a user-chosen folder
% using the multi-line "computeArea" function, which returns an array
% of areas (one per line).
%
%   1) Asks user to select folder containing .fig files
%   2) For each .fig, calls computeArea(...) => areaVals (1..N lines)
%   3) Collects results in a table with columns {Filename, LineIndex, Area}
%   4) Writes table to "AreaResults.xlsx" in that folder
%
% Requirements:
%   - "computeArea" must be on the MATLAB path or in the same folder
%   - "computeArea" can handle multiple lines in each .fig file,
%     returning a vector of areas, one per line.

    % 1) Prompt user for folder
    folderPath = uigetdir(pwd, 'Select Folder Containing .fig Files');
    if isequal(folderPath, 0)
        disp('No folder selected. Exiting.');
        return;
    end

    % 2) Gather all *.fig files in that folder
    figList = dir(fullfile(folderPath, '*.fig'));
    if isempty(figList)
        warning('No .fig files found in folder: %s', folderPath);
        return;
    end

    % We'll collect results in a cell array with 3 columns:
    % {Filename, LineIndex, Area}
    results = {};
    rowCount = 0;

    numFiles = numel(figList);

    % 3) Loop over each .fig
    for iFile = 1:numFiles
        figName = figList(iFile).name;
        figFullPath = fullfile(folderPath, figName);

        % "computeArea" now returns an array of areas, one per line
        % For example, areaVals = [A1, A2, ...]
        % We'll store each line's area on a separate row
        areaVals = computeArea(figFullPath, true);

        % If areaVals is scalar => single line in that fig
        % If vector => multiple lines => multiple rows
        for iLine = 1:numel(areaVals)
            rowCount = rowCount + 1;
            results{rowCount,1} = figName;     % filename
            results{rowCount,2} = iLine;       % line index
            results{rowCount,3} = areaVals(iLine);  % area
        end

        fprintf('Processed %s => found %d lines.\n', figName, numel(areaVals));
    end

    % 4) Convert to table with headings {Filename, LineIndex, Area}
    if rowCount < 1
        warning('No lines/areas found in all .fig files?');
        return;
    end

    T = cell2table(results, ...
        'VariableNames', {'Filename','LineIndex','Area'});

    % Write table => "AreaResults.xlsx"
    outXlsx = fullfile(folderPath, 'AreaResults.xlsx');
    writetable(T, outXlsx);
    fprintf('Wrote results to %s\n', outXlsx);

    disp('Done processing all .fig files.');
end
