function plot_overlapping_cells(varargin)
smooth = load(varargin{1});
num_cells = varargin{2};
save_fname = varargin{3};
t = smooth(1, :);
for i = 1 : num_cells
	c = smooth(1 + i, :);
	plot(t, c, 'color', rand(1, 3));
	hold on 
end
saveas(gcf, save_fname);
exit();
