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