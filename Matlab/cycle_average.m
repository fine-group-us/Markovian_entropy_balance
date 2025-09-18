function [out] = cycle_average(v, Nstep)
% Averages a time series over its cycle: returns an array with Nstep rows.
% v     : input array (N samples × M variables)
% Nstep : cycle length (number of phases, used as modulo)

    dim=size(v);
    N = dim(1);
    out = zeros( Nstep,dim(2));
    
    for i = 1:N
        idx = mod(i-1, Nstep)+1 ; % MATLAB usa índices desde 1
        out(idx,:) = out(idx,:) + v(i,:);
    end

end