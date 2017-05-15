function plot_det_sto_data (varargin)
det = load(varargin{1});
sto = load(varargin{2});
cell_index = varargin{3};
save_file_name = varargin{4};
dt = det(:,1);
dmh1 = det(:,2);
st = sto(1,:);
plot(dt, dmh1, 'b', st, sto(cell_index + 2, :), 'r');
saveas(gcf, save_file_name);
exit;
