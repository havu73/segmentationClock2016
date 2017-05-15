function plot_overlapping_stochastic (varargin)
nargs = length(varargin);
num_round = nargs - 1;
for i = 2:num_round
	raw_file = load(varargin{i});
	t = raw_file(1,:);
	mh1 = raw_file(2,:);
	plot(t, mh1, 'color', rand(1,3));
	hold on
end
saveas(gcf, varargin{nargs});
exit();
