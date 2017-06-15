function plot_smooth_data (varargin)
raw = load (varargin{1});
smooth = load (varargin{2});
cell_index = varargin{3};
pid = fopen (varargin{4});
tid = fopen (varargin{5});
save_file_name = (varargin{6});
% read peak file
p = textscan (pid, '%s', 'delimiter', '\n');
peaks = p {1}{cell_index};
peaks = str2num (peaks);
% read trough file
t = textscan (tid, '%s', 'delimiter', '\n');
troughs = t {1}{cell_index};
troughs = str2num (troughs);
% tranpose the maxtrix of raw and smooth. Now they are, by columns: time, cell 1, cell2, cell3, cell4,...
raw = transpose(raw);
smooth = transpose(smooth);
% time
t_us = raw(:, 1);
t_s = smooth (:, 1);
maxy = max (raw(:, cell_index + 1));
display(maxy);
[m n] = size (peaks);
for i = 1 : n
    p_i = peaks (1, i);
    time = t_s(p_i, 1);
    plot ([time time], [0 maxy], 'm'); 
    legend('hide');
    hold on
end 
[m n] = size (troughs)
for i = 1: n
    t_i = troughs(1, i);
    time = t_s(t_i, 1);
    plot ([time time], [0 maxy], 'g');
    legend('hide');
    hold on
end
hold on
% plot raw data
r = plot (t_us, raw(:, cell_index + 1), 'b')
hold on
% plot smooth data
s = plot (t_s, smooth(:, cell_index + 1), 'r') 
lgd = legend([r s], {'Stochastic Simulation', 'Smoothened Simulation'}, 'FontSize', 14);
xlabel('Time', 'FontSize', 20);
ylabel('mRNA her1 levels', 'FontSize', 20);
saveas(gcf, save_file_name);
exit();
   
% saveas (fig, 'trialsave.png');
% close (fig);
% example call: plot_smooth_data('data01252017/mHer1.txt', 'data01252017/mHer1_smoothed.txt', 1, 'data01252017/mHer1_peaks.txt', 'data01252017/mHer1_troughs.txt');


    
