clear; clc;
theo = processTheoreticalData('theo_5degs.mat', 'theo_5degs_r.mat', true);
% Load experimental data:
load('t8_FT_clipped_avg.mat');   
load('t8_clipped.mat');          

phi = 5;
angle = t8;
force = t8_FT;

force.Fx_r = force.Fx.*cosd(phi) - force.Fz.*sind(phi);
force.Fz_r = force.Fx.*sind(phi) + force.Fz.*cosd(phi);

% theo.Fx_r = theo.Fx.*cos(phi) - theo.Fz.*sin(phi);
% theo.Fz_r = theo.Fx.*sin(phi) + theo.Fz.*cos(phi);
theo.Fx_r = theo.Fx;
theo.Fz_r = theo.Fz;
% row_idx = find(angle.AnteriorLegAngle == 0, 1, 'first');

offset = 180 - angle.AnteriorLegAngle(end);

% % offset = 0;
angle.AnteriorLegAngle = angle.AnteriorLegAngle + offset;


dualAxisFlag = false;
excludeCycles = {'cycle1','cycle14','cycle13','cycle12'};
% excludeCycles = {'cycle1', 'cycle14'};  
% forceOffset.Fx_r = -force.Fx_r(1);
% forceOffset.Fz_r = -force.Fz_r(1);
forceOffset.Fx_r = -0.5;
forceOffset.Fz_r = -0.3;
% forceOffset.Fx_r = 0;
% forceOffset.Fz_r = 0; 

% Generate and save figures:
combineAngleAndForce_saveFigures(angle, force, theo, dualAxisFlag, excludeCycles,forceOffset);


function combineAngleAndForce_saveFigures(angle, force, theo, dualAxisFlag, excludeCycles,forceOffset)

    if nargin < 6
        error('Not enough inputs.');
    end

    %%% 1. Separate angle data into cycles
    cycles = separateAngleCycles(angle);
    cycleNames = fieldnames(cycles);
    if isempty(cycleNames)
        error('No cycle data found. Ensure separateAngleCycles returns cycle data.');
    end
    
    %%% 2. Process each force variable
    forceVars = {'Fx_r', 'Fz_r'};
    for fIdx = 1:numel(forceVars)
        varName = forceVars{fIdx};
        disp(['Processing variable ', varName, '...']);
                
        %%% 2A. Individual Cycle Plots
        for i = 1:length(cycleNames)
            cycName = cycleNames{i};
            cycData = cycles.(cycName);

            % offset = 270 - cycData.AnteriorLegAngle(end);
            % cycData.AnteriorLegAngle = cycData.AnteriorLegAngle - offset;
            
            % 
            % Interpolate force data to cycle time stamps and apply offset only if required.
            % forceOffset.(varName) = -force.(varName)(1);
            yInterp = interp1(force.Time_ms, force.(varName), cycData.Time_ms_, 'spline') + forceOffset.(varName);
                   
            % Combine data and wrap the angle
            combined_data = table(cycData.AnteriorLegAngle, yInterp, 'VariableNames', {'thetaDeg', 'Force'});
            combined_data.thetaDeg_ref = combined_data.thetaDeg + 180;
            idx = (combined_data.thetaDeg_ref > 360);
            combined_data.thetaDeg_ref(idx) = combined_data.thetaDeg_ref(idx) - 360;
            combined_data = sortrows(combined_data, 'thetaDeg_ref');
            
            % ----- Single-Axis Plot for This Cycle -----
            figName = sprintf('%s_%s.fig', cycName, varName);
            figCycle = figure('Name', [varName, ' vs Angle for ', cycName], 'NumberTitle', 'off');
            plot(combined_data.thetaDeg_ref, combined_data.Force, 'o-', 'DisplayName', 'Measured');
            hold on; grid on;
            if ismember(varName, theo.Properties.VariableNames)
                % Wrap theoretical theta if needed (theoretical data is not offset)
                theo_wrapped = wrapTheoretical(theo);
                plot(theo_wrapped.thetaDeg_ref, theo.(varName), 'LineWidth', 1.5, 'DisplayName', 'Theoretical');
            end
            hold off;
            xlabel("Motor Angle (\theta_m)");
            ylabel(varName);
            title([varName, " vs \theta_m for ", cycName]);
            legend('Location', 'bestoutside');
            savefig(figCycle, figName);
            saveas(figCycle, strrep(figName, '.fig', '.jpeg'));
            close(figCycle);
            
            % ----- Dual-Axis Plot for This Cycle (if requested) -----
            if dualAxisFlag
                figDualName = sprintf('DualYAxis_%s_%s.fig', cycName, varName);
                figCycleDual = figure('Name', [varName, ' Dual Y-Axis for ', cycName], 'NumberTitle', 'off');
                
                % Left axis: Experimental data with offset if applicable
                yyaxis left;
                hold on; grid on;
                plot(combined_data.thetaDeg_ref, combined_data.Force, 'o-', 'DisplayName', [cycName, ' (measured)']);
                ylabel(['Experimental ', varName]);
                
                % Right axis: Theoretical data (no offset)
                yyaxis right;
                hold on;
                if ismember(varName, theo.Properties.VariableNames)
                    theo_wrapped = wrapTheoretical(theo);
                    plot(theo_wrapped.thetaDeg_ref, theo.(varName), 'k-', 'LineWidth', 2, 'DisplayName', 'Theoretical');
                end
                ylabel(['Theoretical ', varName]);
                xlabel("Motor Angle (\theta_m)");
                title([varName, ' Dual Y-Axis Plot for ', cycName]);
                legend('Location', 'bestoutside');
                savefig(figCycleDual, figDualName);
                saveas(figCycleDual, strrep(figDualName, '.fig', '.jpeg'));
                close(figCycleDual);
            end
        end
        
        %%% 2B. Overall Overlay Plots (only include cycles not in excludeCycles)
        validCycleNames = setdiff(cycleNames, excludeCycles);
        
        % ----- Overall Single-Axis Overlay Plot -----
        figAllName = sprintf('AllCycles_%s.fig', varName);
        figAll = figure('Name', [varName, ' AllCycles'], 'NumberTitle', 'off');
        hold on; grid on;
        for i = 1:length(validCycleNames)
            cycName = validCycleNames{i};
            cycData = cycles.(cycName);
            
            yInterp = interp1(force.Time_ms, force.(varName), cycData.Time_ms_, 'spline') + forceOffset.(varName);
            
            combined_data = table(cycData.AnteriorLegAngle, yInterp, 'VariableNames', {'thetaDeg', 'Force'});
            combined_data.thetaDeg_ref = combined_data.thetaDeg + 180;
            idx = (combined_data.thetaDeg_ref > 360);
            combined_data.thetaDeg_ref(idx) = combined_data.thetaDeg_ref(idx) - 360;
            combined_data = sortrows(combined_data, 'thetaDeg_ref');
            plot(combined_data.thetaDeg_ref, combined_data.Force, 'o-', 'DisplayName', [cycName, '(measured)']);
        end
        if ismember(varName, theo.Properties.VariableNames)
            theo_wrapped = wrapTheoretical(theo);
            plot(theo_wrapped.thetaDeg_ref, theo.(varName), 'k-', 'LineWidth', 2, 'DisplayName', 'Theoretical');
        end
        hold off;
        xlabel("Motor Angle (\theta_m)");
        ylabel(varName);
        title(['All Cycles + ', varName]);
        legend('Location', 'bestoutside');
        savefig(figAll, figAllName);
        saveas(figAll, strrep(figAllName, '.fig', '.jpeg'));
        close(figAll);
        
        % ----- Overall Dual-Axis Overlay Plot (if requested) -----
        if dualAxisFlag
            figSeparateName = sprintf('SeparateYAxis_%s.fig', varName);
            figSeparate = figure('Name', [varName, ' Separate Y-Axis Plot'], 'NumberTitle', 'off');
            hold on; grid on;
            % Left axis: Experimental data with offset if applicable
            yyaxis left;
            for i = 1:length(validCycleNames)
                cycName = validCycleNames{i};
                cycData = cycles.(cycName);
                
                yInterp = interp1(force.Time_ms, force.(varName), cycData.Time_ms_, 'spline') + forceOffset.(varName);

                combined_data = table(cycData.AnteriorLegAngle, yInterp, 'VariableNames', {'thetaDeg', 'Force'});
                combined_data.thetaDeg_ref = combined_data.thetaDeg + 180;
                idx = (combined_data.thetaDeg_ref > 360);
                combined_data.thetaDeg_ref(idx) = combined_data.thetaDeg_ref(idx) - 360;
                combined_data = sortrows(combined_data, 'thetaDeg_ref');
                plot(combined_data.thetaDeg_ref, combined_data.Force, 'o-', 'DisplayName', [cycName, ' (measured)']);
            end
            ylabel(['Experimental ', varName]);
            
            % Right axis: Theoretical data (no offset)
            yyaxis right;
            if ismember(varName, theo.Properties.VariableNames)
                theo_wrapped = wrapTheoretical(theo);
                plot(theo_wrapped.thetaDeg_ref, theo.(varName), 'k-', 'LineWidth', 2, 'DisplayName', 'Theoretical');
            end
            ylabel(['Theoretical ', varName]);
            xlabel("Motor Angle (\theta_m)");
            title([varName, ' Separate Y-Axis Plot']);
            legend('Location', 'bestoutside');
            savefig(figSeparate, figSeparateName);
            saveas(figSeparate, strrep(figSeparateName, '.fig', '.jpeg'));
            close(figSeparate);
        end
    end
    
    disp('Done creating and saving figures for Fx and Fz.');
end

function theo_wrapped = wrapTheoretical(theo)
% wrapTheoretical Creates a wrapped theta variable from the theoretical data.
%
%   theo_wrapped = wrapTheoretical(theo)
%
%   If the theoretical table does not already have the field thetaDeg_ref, this
%   function creates it by shifting thetaDeg_p by 180° and wrapping any values over 360°.
%
    theo_wrapped = theo;
    if ~ismember('thetaDeg_ref', theo.Properties.VariableNames)
        theo_wrapped.thetaDeg_ref = theo.thetaDeg + 180;
        idxTheo = (theo_wrapped.thetaDeg_ref > 360);
        theo_wrapped.thetaDeg_ref(idxTheo) = theo_wrapped.thetaDeg_ref(idxTheo) - 360;
    end
end

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

function cycles = separateAngleCycles(T)
% separateAngleCycles:
%   Takes a table T with columns:
%       T.Time_ms_          (monotonically increasing time)
%       T.AnteriorLegAngle  (an angle that cycles from ~0 up to >=360, then resets)
%
%   1) Adjusts T.Time_ms_ so that T.Time_ms_(1) = 0 for the entire dataset,
%      i.e. the time is shifted by T.Time_ms_(1).
%   2) Splits T into cycle1, cycle2, cycle3, ... based on angle resets
%      (a new cycle is detected when an entry is less than its predecessor).
%      Only the very first row of the first cycle is reset to 0; subsequent
%      cycles keep their absolute time offsets.
%
%   Returns:
%       cycles - A structure with fields 'cycle1', 'cycle2', etc., each
%                containing the corresponding cycle's table.
%
% Example:
%   cycles = separateAngleCycles(T);

    angle = T.AnteriorLegAngle;
    N = height(T);
    if N < 2
        warning('Not enough data points to separate cycles.');
        cycles = [];
        return
    end

    % Step 1: Re-zero the dataset so that the first row is at 0 ms.
    t0 = T.Time_ms_(1);
    T.Time_ms_ = T.Time_ms_ - t0;  % Now row #1 is 0 ms; subsequent rows keep their offsets

    % Step 2: Identify cycle boundaries.
    % A new cycle starts whenever an angle is less than its previous value.
    cycleStarts = 1;  % First cycle starts at row 1.
    for i = 2:N
        if angle(i) < angle(i-1)
            cycleStarts(end+1) = i; %#ok<AGROW>
        end
    end
    cycleStarts(end+1) = N + 1;  % Mark the end as a boundary.

    % Step 3: Split T into cycles and store in a structure.
    numCycles = length(cycleStarts) - 1;
    cycles = struct();
    for c = 1:numCycles
        iStart = cycleStarts(c);
        iEnd   = cycleStarts(c+1) - 1;
        cycData = T(iStart:iEnd, :);
        
        % Save each cycle as a field in the output structure.
        varName = sprintf('cycle%d', c);
        cycles.(varName) = cycData;
        
        fprintf('Created %s with rows %d..%d, angle from ~%.2f to ~%.2f, time range [%.2f..%.2f] ms.\n',...
                varName, iStart, iEnd, angle(iStart), angle(iEnd), cycData.Time_ms_(1), cycData.Time_ms_(end));
    end
end
