clear;
clc;
%% set working directory
%% read data
Tricho_CH4_prod_yearly_sim_mean = readtable("Tricho_CH4_prod_yearly_sim_mean.csv");
Tricho_CH4_prod_yearly_sim_mean_mat = reshape(Tricho_CH4_prod_yearly_sim_mean.CH4_prod, [180, 360]);
%%
lat_vector = -89.5:1:89.5;
lon_vector = -179.5:1:179.5;
[lon_mat, lat_mat] = meshgrid(lon_vector, lat_vector);

%% CH4 prod mean
[Tricho_CH4_prod_yearly_sim_mean_mat_transform, breaks_new, breaks_lable] = ...
    Transform_for_uneven_colorbar_percentile(Tricho_CH4_prod_yearly_sim_mean_mat*1000, 1);

%%
figure
figure_handle = gcf; % Get current figure handle
set(figure_handle, 'Position', [100, 100, 950, 500]); % [left, bottom, width, height]

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
m_contourf(lon_re, lat_re, Tricho_CH4_prod_yearly_sim_mean_mat_transform(:, ind),50,'LineStyle','none','EdgeColor','none');
%colormap(flipud(hot));
c = colorbar;  % attach colorbar to h
set(c,'ytick',breaks_new,'yticklabel',breaks_lable,...
        'tickdir','out', 'FontName', 'Times New Roman', 'Fontsize',25);
set(get(c,'ylabel'),'String',["{\it Trichodesmium} CH_4 production", "(\mumol N m^{-2} year^{-1})"],...
    'Fontsize',25, 'FontName', 'Times New Roman');
m_coast('patch',[.5 .5 .5]);
m_grid('box','fancy','linestyle',':','gridcolor','black','backcolor',[0.9 0.9 0.9],...
     'FontSize',20, 'FontName', 'Times New Roman');
colormap(flipud(redToWhiteColorbarCustom(10)));
clim([1, length(rmmissing(Tricho_CH4_prod_yearly_sim_mean_mat(:)))]);
text(-3, 1.6, 'a', 'FontSize',35, 'FontName', 'Times New Roman')
%%
% print(gcf, 'Figures/fig_CH4_prod_sim_mean', '-dpng', '-r600');
%% Tricho CH4 contribution

Tricho_CH4_contribution_sim_mean = readtable("Tricho_CH4_contribution_sim_mean.csv");
Tricho_CH4_contribution_sim_mean_mat = reshape(Tricho_CH4_contribution_sim_mean.CH4_contri, [180, 360]);

%% Transform_for_uneven_colorbar_percentile
[Tricho_CH4_contribution_sim_mean_mat_transform, breaks_new, breaks_lable] = ...
    Transform_for_uneven_colorbar_percentile(Tricho_CH4_contribution_sim_mean_mat, 1);
%%
figure
figure_handle = gcf; % Get current figure handle
set(figure_handle, 'Position', [100, 100, 950, 500]); % [left, bottom, width, height]

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
m_contourf(lon_re, lat_re, Tricho_CH4_contribution_sim_mean_mat_transform(:, ind),10,'LineStyle','none','EdgeColor','none');
%colormap(cmocean('balance', 10));
%clim([0.8 2]);
colormap(sky(10));
c = colorbar;  % attach colorbar to h
set(c,'ytick',breaks_new,'yticklabel',breaks_lable,...
        'tickdir','out','FontSize',25, 'FontName', 'Times New Roman');
set(get(c,'ylabel'),'String',["Percentage of {\it Trichodesmium}", "contribution to CH_4 flux"],...
    'FontSize',25, 'FontName', 'Times New Roman');
m_coast('patch',[.5 .5 .5]);
m_grid('tickdir','out','linewi',2, 'backcolor', [0.9 0.9 0.9],...
    'FontSize',20, 'FontName', 'Times New Roman');
%text(-3, 1.6, 'd', 'fontsize',25)
%colormap(flipud(m_colmap('green',10)))
clim([1, length(rmmissing(Tricho_CH4_contribution_sim_mean_mat(:)))]);
text(-3, 1.6, 'b', 'fontsize',35, 'FontName', 'Times New Roman')
%%
%print(gcf, 'Figures/fig_CH4_contri_sim_mean', '-dpng', '-r600');