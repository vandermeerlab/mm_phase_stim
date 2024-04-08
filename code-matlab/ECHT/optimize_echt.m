function [optimal_parameters, ga_output] = optimize_echt(csc, filt_phase, fband, Fs, nSamples, popsz, bounds_window)
% optimize_echt
%   optimal_parameters = optimize_echt(csc, filt_phase, fbands, Fs, nSamples, popsz, bounds_window)
%
%     Input:
%         wsz: window size in tenths of a second (for faster optimization)
%         csc: mvdmlab TSD struct
%         filt_phase: array that containns hilbert transformed phases of 'csc'
%         fband: Frequency band with [lfq hfq]
%         Fs: sampling frequency of csc
%         nSamples: Number of samples to be test this method on
%         popsz: population size for the genetic algorithm
%         bounds_window: Upper and  lower limit for window_size in tenths of seconds
%     Output: 
%         optimal_parameters : parameters for the winning run


assert(~isempty(which('ga')), 'genetic algorithm function ga.m not found, is the Global Optimization Toolbox installed?')

problem.options = optimoptions('ga', 'PlotFcn', @gaplotbestf, 'PopulationSize', popsz,'MaxStallGenerations',500,'MaxGenerations',1000, 'UseParallel', true);
problem.solver = 'ga';

ang_var_of_diff = @(x) 1-abs(mean(exp(1i*x{1}')./exp(1i*x{2}')));

% Changing the above to be the mean variance across all freq. bands
% ang_var_of_diff = @(x) mean(1-abs(mean(exp(1i*x{1}')./exp(1i*x{2}'))));

problem.fitnessfcn = @(x) ang_var_of_diff(wrapper_echt(x, csc, filt_phase, fband, Fs, nSamples));
% x(1) = window_length
% x(2) = filter_order
% x(3) = edge
% x(4) = ar_order

problem.nvars = 1;
problem.intcon = 1;
problem.lb = [bounds_window(1)];
problem.ub = [bounds_window(2)];

problem.x0 = ceil(mean([problem.lb; problem.ub]));

% A*x ≤ b
% problem.Aineq = [-1 3 0 0]; %window_length x(1) > 3 * filter_order x(2)
% problem.bineq = -1;

% check if the function evaluates
feval(problem.fitnessfcn, problem.x0);

% run solver
[x,fval,exitflag,ga_output] = ga(problem);

optimal_parameters = [];
optimal_parameters.window_length = x(1);
optimal_parameters.fval = fval;
    
end