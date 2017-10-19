function [TTLPulses] = countTTLPulses (Report)
% [TTLPulses] = countTTLPulses (Report)
% 
% Calculates the unique TTL pulses received by the LiveTrack device. It is
% useful to check if the number corresponds to the expected TR count.
% 
% Input
% Report - struct with data as it was tracked by the liveTrack, in the
% standard format returned by the liveTrack itself.
% 
% Outputs
% TTLPulses - number of the recorded TTL pulses

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
