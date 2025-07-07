clear;
clc;
%% set working directory

%% read Tricho nfix yearly data
% `NA` in csv file exported from R should be first replaced with `NaN` 

mlp_Tricho_nfix_prediction_yearly_summary = ...%  from Prediction_and_Uncertainty_analysis.Rmd 
    readtable("mlp_Tricho_nfix_prediction_yearly_summary.csv"); 
Tricho_nfix_pred_yearly_mean_mat = ...
    reshape(mlp_Tricho_nfix_prediction_yearly_summary.mean, [180, 360]);
Tricho_nfix_pred_yearly_median_mat = ...
    reshape(mlp_Tricho_nfix_prediction_yearly_summary.median, [180, 360]);
%%
lat_vector = -89.5:1:89.5;
lon_vector = -179.5:1:179.5;
[lon_mat, lat_mat] = meshgrid(lon_vector, lat_vector);
%% figure Tricho nfix yearly mean 
[Tricho_nfix_pred_yearly_mean_mat_transform, breaks_new, breaks_lable] = ...
    Transform_for_uneven_colorbar_percentile(Tricho_nfix_pred_yearly_mean_mat, 1);
%%
figure
figure_handle = gcf; % Get current figure handle
set(figure_handle, 'Position', [100, 100, 800, 500]); % [left, bottom, width, height]

% Tighten axes to minimize white space
tightInset = get(gca, 'TightInset');
set(gca, 'Position', [tightInset(1)+0.04, tightInset(2)+0.15, ...
    1-0.15, 1-0.3]);
ind = [211:360 1:210];
lon_re = lon_mat(:, ind);
lon_re(lon_re<0) = lon_re(lon_re<0)+360;
lat_re = lat_mat(:, ind);
lon_re(lon_re > 30) = lon_re(lon_re > 30) - 360;
hold on
m_proj('robinson','lon',[-330 30]); % Specify map domain
m_contourf(lon_re, lat_re, Tricho_nfix_pred_yearly_mean_mat_transform(:, ind),50,'LineStyle','none','EdgeColor','none');
%colormap(flipud(hot));
c = colorbar;  % attach colorbar to h
set(c,'ytick',breaks_new,'yticklabel',breaks_lable,...
        'tickdir','out');
set(get(c,'ylabel'),'String','N_2 fixation (mmol N m^{-2} year^{-1})',...
    'fontsize',20, 'Position', [2.3 max(breaks_new)*0.5]);
m_coast('patch',[.5 .5 .5]);
m_grid('box','fancy','linestyle',':','gridcolor','black','backcolor',[0.9 0.9 0.9],...
     'fontsize',10);
colormap(flipud(m_colmap('green',10)))
clim([0, length(rmmissing(Tricho_nfix_pred_yearly_mean_mat(:)))]);