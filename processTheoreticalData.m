
function theo = processTheoreticalData(inputFile, outputFile, plotFlag)
% processTheoreticalData Processes and saves theoretical data.
%
%   theo = processTheoreticalData(inputFile, outputFile, plotFlag)
%
%   This function loads theoretical data from the specified input file,
%   creates a new reference variable (thetaDeg_ref) by shifting thetaDeg_p by 180°,
%   removes the first (repeated) row, adjusts any values over 360°, sorts the table,
%   and saves the processed data to the output file. Optionally, it generates and saves
%   plots for Fz and Fx.
%
%   Inputs:
%       inputFile  - String, name of the .mat file containing the theoretical data.
%       outputFile - String, name of the .mat file to save the processed data.
%       plotFlag   - Logical, if true (1) the function generates verification plots.
%
%   Output:
%       theo       - Processed theoretical data (table).
%
%   Example:
%       theo = processTheoreticalData('theo_0degs.mat', 'theo_0degs_r.mat', true);
%
% Author: Your Name
% Date:   YYYY-MM-DD

    % Load the input data (assumes the file contains one variable)
    data = load(inputFile);
    fieldNames = fieldnames(data);
    % Assume the first variable in the file contains the theoretical data
    theo = data.(fieldNames{1});
    
    % Create new reference variable: shift thetaDeg_p by 180°
    theo.thetaDeg_ref = theo.thetaDeg_p + 180;
    
    % Remove first repeated element
    theo(1,:) = [];
    
    % Adjust values greater than 360°
    idx = (theo.thetaDeg_ref > 360);
    theo.thetaDeg_ref(idx) = theo.thetaDeg_ref(idx) - 360;
    
    % Sort the table based on the reference variable
    theo = sortrows(theo, 'thetaDeg_ref');
    
    % Save the processed theoretical data to the output file
    save(outputFile, 'theo');
    
    % Optionally, plot the processed data to verify results
    if plotFlag
        % Determine a base name from outputFile for plot files
        [~, baseName, ~] = fileparts(outputFile);
        
        % Plot Fz versus thetaDeg_ref
        fxFig = figure;
        plot(theo.thetaDeg_ref, theo.Fz);
        xlabel("\theta_{ref}");
        ylabel("F_z");
        title('Processed Theoretical Data: F_z');
        savefig(fxFig, [baseName '_Fz.fig']);
        saveas(fxFig, [baseName '_Fz.png']);
        close(fxFig);
        
        % Plot Fx versus thetaDeg_ref
        fzFig = figure;
        plot(theo.thetaDeg_ref, theo.Fx);
        xlabel("\theta_{ref}");
        ylabel("F_x");
        title('Processed Theoretical Data: F_x');
        savefig(fzFig, [baseName '_Fx.fig']);
        saveas(fzFig, [baseName '_Fx.png']);
        close(fzFig);
    end
end
