function [optimal_parameters, ga_output] = optimize_sspe(csc, fbands, filt_phase, Fs, nSamples, popsz, bounds_window, bounds_exp)
% optimize_rez
%   optimal_parameters = optimize_sspe(csc, fbands, filt_phase, Fs, nSamples, popsz, bounds_window, bounds_exp)
%   
%
%     Input:
%         csc: mvdmlab TSD struct
%         fbands: nx2 array of frequency bands [[lfq1 hfq];[lfq1 hfq2]; ..., lfqn hfqn]]
%         filt_phase: cell array where each cell contians hilbert transformed phases of 'csc'
%         Fs: sampling frequency of csc
%         nSamples: Number of samples to be test this method on
%         popsz: population size for the genetic algorithm

%   
%     Output: 
%         optimal_parameters : parameters for the winning run


assert(~isempty(which('ga')), 'genetic algorithm function ga.m not found, is the Global Optimization Toolbox installed?')

problem.options = optimoptions('ga', 'PlotFcn', @gaplotbestf, 'PopulationSize', popsz, 'UseParallel', true);
problem.solver = 'ga';

% ang_var_of_diff = @(x, y) 1-abs(mean(exp(1i*x)./exp(1i*y)));

% Changing the above to be the mean variance across all freq. bands
ang_var_of_diff = @(x) mean(1-abs(mean(exp(1i*x{1}')./exp(1i*x{2}'))));

% wrapper_sspe(wsz, fq_exp, csc, fbands, filt_phase, Fs, nSamples)
problem.fitnessfcn = @(x) ang_var_of_diff(wrapper_sspe(x(1), [x(2) x(3) x(4) x(5)], csc, fbands, filt_phase, Fs, nSamples));
% x(1) = window_length
% x(2) = fq(1)
% x(3) = fq(2)
% x(4) = fq(3)
% x(5) = fq(4)

problem.nvars = 5;
problem.intcon = 1;
problem.lb = [bounds_window(1) 10^(-bounds_exp(2)) 10^(-bounds_exp(2)) 10^(-bounds_exp(2)) 10^(-bounds_exp(2))];
problem.ub = [bounds_window(2) 10^(-bounds_exp(1)) 10^(-bounds_exp(1)) 10^(-bounds_exp(1)) 10^(-bounds_exp(1))];

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
optimal_parameters.fq_exp = [x(2) x(3) x(4) x(5)];
optimal_parameters.fval = fval;
    
end