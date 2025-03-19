
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
