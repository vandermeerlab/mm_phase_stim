function [angle_tsd, wheel_tsd, bad_jumps] = ConvertQEStatesToAngle(cfg_in, state_tsd)


cfg_def = [];
cfg_def.dphi = 1; % encoder angle step -- property of the encoder (should be obtained by calibration procedure)
cfg_def.angle_to_turn  = 1024; % this is the conversion factor for coverting integer "angle" values into # of wheel turns. This was determined empirically with a recording session on 3/3/2022. 


% state definitions:
% high-high: 1
% high-low: 2
% low-high: 3
% low-low: 4

% this variable maps state transitions to forward (1) or reverse (-1) steps
cfg_def.tt = [0 -1 1 NaN; 1 0 NaN -1; -1 NaN 0 1; NaN 1 -1 0];
% cfg_def.tt = [0 -1 1 0; 1 0 0 -1; -1 0 0 1; 0 1 -1 0];
% example: given a transition from high-high (1) to high-low (2),
% cfg_def.tt(1, 2) is -1, so reverse step
  
cfg = ProcessConfig(cfg_def, cfg_in);

angle_tsd = state_tsd;
angle_tsd.data = NaN*angle_tsd.data;
angle_tsd.data(1) = 0; % assume we start at 0


% loop over state inputs; look up step given transition between current
% state and previous state at each step
for iState = 2:length(state_tsd.data)   
    angle_tsd.data(iState) = cfg.tt(state_tsd.data(iState-1), state_tsd.data(iState));
    
end    % units of "angle" here are actual integer values that represent identical forward or backward increments of the wheel (i.e. not an actual angle value)

% angle_tsd.data can contain NaNs in case of illegal transitions
bad_jumps = isnan(angle_tsd.data);

% Linearly interpolate these transitions
angle_tsd.data(bad_jumps) = interp1(angle_tsd.tvec(~bad_jumps), angle_tsd.data(~bad_jumps), angle_tsd.tvec(bad_jumps));
% gstart = 1; gend = length(angle_tsd.data);
% if ~isempty(bad_jumps)
%     for iS = 1:length(bad_jumps)
%         if gstart > length(angle_tsd.data)
%             break;
%         end
%         gend = min(gend, bad_jumps(iS)) - 1;
%         % Sanity check
%         assert(sum(isnan(angle_tsd.data(gstart:gend)))==0);
%         angle_tsd.data(gstart:gend) = cumsum(angle_tsd.data(gstart:gend));
%         gstart = bad_jumps(iS)+1;
%     end
% else    
%     
% end


angle_tsd.data = cumsum(angle_tsd.data);
wheel_tsd.tvec = angle_tsd.tvec;  
wheel_tsd.data = angle_tsd.data./cfg.angle_to_turn; 