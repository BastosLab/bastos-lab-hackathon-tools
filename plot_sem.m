function [] = plot_sem(xs, ys, sems, color, trial_axis)
%PLOT_SEM Summary of this function goes here
%   Detailed explanation goes here
if ~exist('sems', 'var')
    sems = 2;
end
if ~exist('color', 'var')
    color = 'k';
end
if ~exist('trial_axis', 'var')
    trial_axis = size(size(ys));
    trial_axis = trial_axis(end);
end
ymean = squeeze(nanmean(ys, trial_axis));
ysem = squeeze(nanstd(ys,[],trial_axis) / sqrt(size(ys, trial_axis)));
fill([xs; flipud(xs)], [ymean - sems * ysem; flipud(ymean + sems * ysem)], ...
    color, 'FaceAlpha', 0.25, 'linestyle', 'none');
line(xs,ymean, 'Color', color);
end

