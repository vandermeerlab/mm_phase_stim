function [optimal_parameters ga_output] = optimize_echt(csc, nSamples, fbands, Fs, seed, bounds_window)
% optimize_echt
%   optimal_parameters = phastimate_optimize(epochs, truephase, filter_objects)
%
%   Input:
%     data is a time x epoch matrix
%     truephase is the true phase in radians at the end of the epoch
%     filter_objects_by_order is a cell array of digitalFilter objects indexed by order
%     bounds_NN is lower and upper bound, e.g. [160 500]
%
%   Out: 
%     optimal_parameters : parameters for the winning run

%TODO: make population_size a parameter


assert(~isempty(which('ga')), 'genetic algorithm function ga.m not found, is the Global Optimization Toolbox installed?')

problem.options = optimoptions('ga', 'PlotFcn', @gaplotbestf, 'PopulationSize', 100);
problem.solver = 'ga';

% [true_phase, output_phase] = wrapper_echt(csc, wsz, fbands, Fs, nSamples, seed)


% ang_var_of_diff = @(x, y) 1-abs(mean(exp(1i*x)./exp(1i*y)));

% Changing the above to be the mean variance across all freq. bands
ang_var_of_diff = @(x) mean(1-abs(mean(exp(1i*x{1}')./exp(1i*x{2}'))));

problem.fitnessfcn = @(x) ang_var_of_diff(wrapper_echt(csc, x, fbands, Fs, nSamples, seed));
% x(1) = window_length
% x(2) = filter_order
% x(3) = edge
% x(4) = ar_order

problem.nvars = 1;
% problem.intcon = 1:4;
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