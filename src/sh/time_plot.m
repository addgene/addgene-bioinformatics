
t_spades = load('time_A11967A_sW0154_FASTQ_spades.out'); % [s]
t_apc    = load('time_A11967A_sW0154_FASTQ_apc.out');    % [s]
s_inp    = load('time_A11967A_sW0154_FASTQ_size_inp.out');   % [MB]
s_out    = load('time_A11967A_sW0154_FASTQ_size_out.out');   % [MB]

[n_total, x_total] = hist(t_spades + t_apc, 20);

figure(1);
bar(x_total, n_total);

limits = axis;
x1 = limits(1); x2 = limits(2); xW = x2 - x1;
y1 = limits(3); y2 = limits(4); yW = y2 - y1;

gct = [];
gct = [gct, text(x2 - xW * 0.45, y2 - yW * 0.05, sprintf('Total time: %.1f hrs', sum(t_spades + t_apc) / 60.0 / 60.0))];
gct = [gct, text(x2 - xW * 0.45, y2 - yW * 0.10, sprintf('Mean time: %.1f min', mean(t_spades + t_apc) / 60))];
gct = [gct, text(x2 - xW * 0.45, y2 - yW * 0.15, sprintf('Input disk usage: %d MB', s_inp))];
gct = [gct, text(x2 - xW * 0.45, y2 - yW * 0.20, sprintf('Output disk usage: %d MB', s_out))];
set(gct, 'FontName', 'Helvetica', 'FontSize', 14, 'FontWeight', 'demi');

gcl = [];
gcl = [gcl, xlabel('Total processing time [s]')];
gcl = [gcl, ylabel('Number of FASTQ file pairs')];
set(gcl, 'FontName', 'Helvetica', 'FontSize', 16, 'FontWeight', 'bold');

set(gca, 'FontName', 'Helvetica', 'FontSize', 16, 'FontWeight', 'bold');

print('-dpng', 'time_plot.png');
