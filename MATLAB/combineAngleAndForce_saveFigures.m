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
