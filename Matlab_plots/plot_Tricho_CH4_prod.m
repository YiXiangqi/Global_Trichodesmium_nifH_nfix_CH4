clear;
clc;
%% set working directory

%% read data
Tricho_CH4_prod_yearly_mean = readtable("Tricho_CH4_prod_yearly_mean.csv");
Tricho_CH4_prod_yearly_mean_mat = reshape(Tricho_CH4_prod_yearly_mean.CH4_prod, [180, 360]);
%%
lat_vector = -89.5:1:89.5;
lon_vector = -179.5:1:179.5;
[lon_mat, lat_mat] = meshgrid(lon_vector, lat_vector);

%% CH4 prod mean
[Tricho_CH4_prod_yearly_mean_mat_transform, breaks_new, breaks_lable] = ...
    Transform_for_uneven_colorbar_percentile(Tricho_CH4_prod_yearly_mean_mat*1000, 1);

%%
% Set figure size
figure
figure_handle = gcf; % Get current figure handle
set(figure_handle, 'Position', [100, 100, 600, 350]); % [left, bottom, width, height]

% Tighten axes to minimize white space
tightInset = get(gca, 'TightInset');
set(gca, 'Position', [tightInset(1)+0.02, tightInset(2)+0.15, ...
    1-0.15, 1-0.35]);
ind = [211:360 1:210];
lon_re = lon_mat(:, ind);
lon_re(lon_re<0) = lon_re(lon_re<0)+360;
lat_re = lat_mat(:, ind);
lon_re(lon_re > 30) = lon_re(lon_re > 30) - 360;
hold on
m_proj('robinson','lon',[-330 30]); % Specify map domain
m_contourf(lon_re, lat_re, Tricho_CH4_prod_yearly_mean_mat_transform(:, ind),50,'LineStyle','none','EdgeColor','none');
%colormap(flipud(hot));
c = colorbar;  % attach colorbar to h
set(c,'ytick',breaks_new,'yticklabel',breaks_lable,...
        'tickdir','out', 'FontName','Times New Roman');
set(get(c,'ylabel'),'String',["{\it Trichodesmium} CH_4 production", "(\mumol N m^{-2} year^{-1})"],...
    'fontsize', 14,'FontName','Times New Roman');

m_coast('patch',[.5 .5 .5]);
m_grid('box','fancy','linestyle',':','gridcolor','black','backcolor',[0.9 0.9 0.9],...
     'fontsize',12, 'FontName','Times New Roman');
colormap(flipud(redToWhiteColorbarCustom(10)));
%colormap(flipud(magma(10)));
clim([1, length(rmmissing(Tricho_CH4_prod_yearly_mean_mat(:)))]);
text(-3.3, 1.6, 'b', 'FontSize', 20, 'FontWeight','bold', 'FontName','Times New Roman');
%%
%print(gcf, 'Figures/fig_3b', '-dpng', '-r600');
%% Global CH4 flux
%% CH4 flux mapping
% reading data
% ocean_ch4.nc is from Weber et al.(2019, Nature communications)
CH4_flux = transpose(ncread("ocean_ch4.nc", "Fch4_diffusive"));
CH4_flux_lon = ncread("ocean_ch4.nc", "LON");
CH4_flux_lon(CH4_flux_lon>180) = CH4_flux_lon(CH4_flux_lon>180) - 360;
CH4_flux_lat = ncread("ocean_ch4.nc", "LAT");
[CH4_flux_lon_mat, CH4_flux_lat_mat] = meshgrid(CH4_flux_lon, CH4_flux_lat);
CH4_flux_lon_mat = CH4_flux_lon_mat(:, [721:1440 1:720]);
CH4_flux = CH4_flux(:, [721:1440 1:720]);
%%
lat_vector = -89.5:1:89.5;
lon_vector = -179.5:1:179.5;
[lon_mat, lat_mat] = meshgrid(lon_vector, lat_vector);
CH4_flux_regrid = interp2(CH4_flux_lon_mat, CH4_flux_lat_mat, CH4_flux, lon_mat, lat_mat);
%writematrix(CH4_flux_regrid(:), "Global_CH4_flux.csv") % for analysis in R server
%% figure global CH4 flux
[CH4_flux_regrid_transform, breaks_new, breaks_lable] = ...
    Transform_for_uneven_colorbar_percentile(CH4_flux_regrid, 2);
%%
figure
ind = [211:360 1:210];
lon_re = lon_mat(:, ind);
lon_re(lon_re<0) = lon_re(lon_re<0)+360;
lat_re = lat_mat(:, ind);
lon_re(lon_re > 30) = lon_re(lon_re > 30) - 360;
hold on
m_proj('robinson','lon',[-330 30]); % Specify map domain
m_contourf(lon_re, lat_re, CH4_flux_regrid_transform(:, ind),10,'LineStyle','none','EdgeColor','none');
colormap(cmocean('balance', 11));
%clim([0.8 2]);
%colormap(flipud(hot));
c = colorbar;  % attach colorbar to h
set(c,'ytick',breaks_new,'yticklabel',breaks_lable,...
        'tickdir','out');
set(get(c,'ylabel'),'String','CH4 flux',...
    'fontsize',20, 'Position', [2.3 max(breaks_new)*0.5]);
m_coast('patch',[.5 .5 .5]);
m_grid('tickdir','out','linewi',2, 'backcolor', [0.95 0.95 0.95]);
%text(-3, 1.6, 'd', 'fontsize',25)
%colormap(flipud(m_colmap('green',10)))
clim([0, length(rmmissing(CH4_flux_regrid(:)))]);

%% Tricho CH4 contribution

Tricho_CH4_contri_yearly_mean = readtable("Tricho_CH4_contri_yearly_mean.csv");
Tricho_CH4_contri_yearly_mean_mat = reshape(Tricho_CH4_contri_yearly_mean.Tricho_contri, [180, 360]);

%% Transform_for_uneven_colorbar_percentile
[Tricho_CH4_contri_yearly_mean_mat_transform, breaks_new, breaks_lable] = ...
    Transform_for_uneven_colorbar_percentile(Tricho_CH4_contri_yearly_mean_mat, 1);
%%
figure
figure_handle = gcf; % Get current figure handle
set(figure_handle, 'Position', [100, 100, 600, 350]); % [left, bottom, width, height]

% Tighten axes to minimize white space
tightInset = get(gca, 'TightInset');
set(gca, 'Position', [tightInset(1)+0.02, tightInset(2)+0.15, ...
    1-0.15, 1-0.35]);
ind = [211:360 1:210];
lon_re = lon_mat(:, ind);
lon_re(lon_re<0) = lon_re(lon_re<0)+360;
lat_re = lat_mat(:, ind);
lon_re(lon_re > 30) = lon_re(lon_re > 30) - 360;
hold on
m_proj('robinson','lon',[-330 30]); % Specify map domain
m_contourf(lon_re, lat_re, Tricho_CH4_contri_yearly_mean_mat_transform(:, ind),10,'LineStyle','none','EdgeColor','none');
colormap(sky(10))
%colormap(flipud(sky(10)))
%colormap(cmocean('balance', 10));
%clim([0.8 2]);
%colormap(flipud(hot));
c = colorbar;  % attach colorbar to h
set(c,'ytick',breaks_new,'yticklabel',breaks_lable,...
        'tickdir','out', 'FontName','Times New Roman');
set(get(c,'ylabel'),'String',["Percentage of {\it Trichodesmium}", "contribution to CH_4 flux"],...
    'fontsize',14, 'FontName','Times New Roman');
m_coast('patch',[.5 .5 .5]);
m_grid('tickdir','out','linewi',2, 'backcolor', [0.9 0.9 0.9],...
    'fontsize',12, 'FontName','Times New Roman');
%text(-3, 1.6, 'd', 'fontsize',25)
%colormap(flipud(m_colmap('green',10)))
clim([0, length(rmmissing(Tricho_CH4_contri_yearly_mean_mat(:)))]);
text(-3.3, 1.6, 'c', 'FontSize', 20, 'FontWeight','bold', 'FontName','Times New Roman');
%%
%print(gcf, 'Figures/fig_3c', '-dpng', '-r600');