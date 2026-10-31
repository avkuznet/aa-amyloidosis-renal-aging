function s = funflare(t,Slo,Shi,period,duty,ton,w)
% Smoothed periodic inflammatory flare profile for total SAA (uM).
% Transitions are smoothed over a width w (s), reflecting the 1-2 day
% rise and fall times of the acute-phase SAA response, which also keeps the
% right-hand side differentiable for the adaptive solver.
ph = mod(t,period);
sig = @(x) 1./(1+exp(-x));
pulse = sig((ph-ton)/w) - sig((ph-ton-duty*period)/w);
s = Slo + (Shi-Slo)*pulse;
end
