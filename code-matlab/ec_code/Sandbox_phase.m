hold on
n_phases = 5;
c_ord = linspecer(5);
x = -pi:pi/(50*n_phases):pi;
wave_phase = sin(x);
for iPhase = 1:n_phases
    
    plot(-99+(100*iPhase):(100*iPhase), wave_phase(-99+(100*iPhase):(100*iPhase)), 'color', c_ord(iPhase,:), 'linewidth', 5)
    
end
% axis off
%         title(sprintf('%.1f to %.1f Hz', f_list{iF}(1), f_list{iF}(2)), 'fontsize', font_size)

%% try an overlapping 

polar(rad2deg(wave_phase))


%%
        n_phases = 8;
        fprintf('\nPhase split into %1d bins\n', n_phases)
        [~, base_edges, ~] = histcounts(-pi:pi, n_phases, 'BinLimits', [-pi, pi]);
        
        
        % circ_shift version
%         edges_shift = circshift(base_edges, [2,1]);
        edges_shift = base_edges(2:2:end);
        
        edges = base_edges(1:2:end-1);
        
        
        % same thing but shifted
                n_phases = 4;
        fprintf('\nPhase split into %1d bins\n', n_phases)
        [~, edges_shift, ~] = histcounts(-pi +(pi/4):(pi+pi/4), n_phases, 'BinLimits', [(-pi +(pi/4)), (pi+pi/4)]);