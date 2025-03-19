function readCsv()
    
    % Step 1: Let user pick a .csv file
    [filename, pathname] = uigetfile('*.csv','Select a CSV file');
    
    if isequal(filename,0)
        disp('User canceled file selection. No file imported.');
        return;
    end

    % Construct full path
    csvFile = fullfile(pathname, filename);

    % Step 2: Read the CSV, first row as column headers
    % (Adapt 'Delimiter' if needed for other separators)
    T = readtable(csvFile, 'ReadVariableNames', true);

    % Step 3: Preview the first few rows
    disp('Preview of imported data:');
    head(T);

    % Step 4: Ask user for .mat filename
    outMat = input('Enter output MAT filename (e.g. "myData.mat"): ','s');
    if isempty(outMat)
        outMat = 'ImportedData.mat';  % default name if user leaves blank
    end
    
    
    % Step 5: Save the table T in that .mat file
    save(outMat, 'T');
    fprintf('Table variable "T" saved to %s\n', outMat);
end
