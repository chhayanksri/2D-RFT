function areaVals = computeArea(figFile, saveFlag)
% computeArea  Computes signed area(s) under one or more x-y curves from a .fig file,
% where y<0 region is negative area (red), y>=0 is positive area (green).
%
% If multiple line objects exist, it calculates an area for each line
% and returns them in an array areaVals. A single new figure is produced
% showing shading for all lines together.
%
% areaVals = computeArea(figFile, saveFlag)
%
% Example:
%   A = computeArea('myPlot.fig', true);
%   % A is [A1; A2; ...] for each line. The new figure is saved if saveFlag = true.
%
% Steps:
%   1) Open old figure (invisible).
%   2) Find axes, lines. Extract old figure's title string.
%   3) For each line: 
%       (a) sort by x,
%       (b) split by sign of y => y<0 vs y>=0,
%       (c) integrate each portion => areaNeg + areaPos => net area
%   4) Create one new figure. For each line, shade negative portion in red,
%      positive portion in green, and plot the line in black.
%   5) Title is set to the original figure's title. 
%   6) If saveFlag, save new figure as e.g. "myPlot_shaded.fig"/.png
%   7) areaVals is returned as a column vector, one area per line.
%

    if nargin < 2
        saveFlag = false;  % default
    end

    % 1) Open old figure invisibly
    oldFigHandle = openfig(figFile, 'invisible');

    % 2) Find axes & line objects
    axHandles = findobj(oldFigHandle, 'Type', 'axes');
    if isempty(axHandles)
        error('No axes found in the figure.');
    elseif numel(axHandles) > 1
        warning('Multiple axes found; using the first one.');
    end
    ax = axHandles(1);

    % Extract old figure's title string
    oldTitleStr = ax.Title.String; 

    % Find all lines in that axes
    linesInFig = findobj(ax, 'Type','line');
    if isempty(linesInFig)
        error('No line objects found in the figure.');
    end

    numLines = numel(linesInFig);
    fprintf('Found %d line objects.\n', numLines);

    % We'll store each line's area in areaVals
    areaVals = zeros(numLines,1);

    % 3) Close the original figure (we won't use it further)
    % close(oldFigHandle);

    % 4) Create new figure for shading all lines
    figShade = figure('Name','Signed Area by Y-Sign','NumberTitle','off');
    hold on; grid on;

    for iLine = 1:numLines
        ln = linesInFig(iLine);

        % Extract & sort
        x = ln.XData;
        y = ln.YData;
        [x, idx] = sort(x);
        y = y(idx);

        % Split by sign of y
        idxNeg = (y<0);
        idxPos = (y>=0);

        xNeg = x(idxNeg);
        yNeg = y(idxNeg);

        xPos = x(idxPos);
        yPos = y(idxPos);

        % Integrate
        areaNeg = trapz(xNeg,yNeg);  % y<0 portion => typically negative
        areaPos = trapz(xPos,yPos);  % y>=0 portion => typically positive
        areaVal = areaPos + areaNeg; % effectively subtract negative portion

        areaVals(iLine) = areaVal;

        % 5) Shading for this line:
        % Negative portion => red
        if ~isempty(xNeg)
            patchX_neg = [xNeg, fliplr(xNeg)];
            patchY_neg = [yNeg, fliplr(zeros(size(yNeg)))];
            patch(patchX_neg, patchY_neg,...
                  'r','FaceAlpha',0.3,'EdgeColor','none');
        end

        % Positive portion => green
        if ~isempty(xPos)
            patchX_pos = [xPos, fliplr(xPos)];
            patchY_pos = [zeros(size(yPos)), fliplr(yPos)];
            patch(patchX_pos, patchY_pos,...
                  'g','FaceAlpha',0.3,'EdgeColor','none');
        end

        % Plot the curve in black
        plot(x, y, 'k-','LineWidth',1.5);

        fprintf('Line #%d: areaVal = %g\n', iLine, areaVal);
    end

    xlabel('X');
    ylabel('Y');
    title(oldTitleStr);
    legend off;  % or keep off if there's too many patches
    axis tight;

    % 6) If saveFlag, save new figure
    if saveFlag
        [~, baseName, ~] = fileparts(figFile);
        figOutName = sprintf('%s_shaded.fig', baseName);
        pngOutName = sprintf('%s_shaded.png', baseName);
        savefig(figShade, figOutName);
        saveas(figShade, pngOutName);
        fprintf('Saved shading figure as "%s" and "%s"\n', figOutName, pngOutName);
    end
end
