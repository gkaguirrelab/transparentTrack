function [TTLPulses] = CountTTLPulses (Report)
% Calculates the unique TTL pulses received by the LiveTrack device. It is
% useful to check if the number corresponds to the TR count.

allPulses = find ([Report.Digital_IO1]);
if isempty(allPulses)
    TTLPulses = 0;
else
    spacing = diff(allPulses);
    for ii = 1:length(spacing)
        if spacing (ii) == 1
            adjacent (ii) = 1;
        else
            adjacent (ii) = 0;
        end
    end
    TTLPulses = length (allPulses) - length (find(adjacent));
end
