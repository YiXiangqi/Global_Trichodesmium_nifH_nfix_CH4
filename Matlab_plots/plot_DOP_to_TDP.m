clear;
clc;
%% set working directory

%% WOA 2013  SRP mapping
% reading data
SRP_path_WOA13 = 'woa13_all_p00_01.nc';
SRP_WOA13_lon = ncread(SRP_path_WOA13,"lon");
SRP_WOA13_lon = repmat(SRP_WOA13_lon', 180, 1);
SRP_WOA13_lat = ncread(SRP_path_WOA13,"lat");
SRP_WOA13_lat = repmat(SRP_WOA13_lat, 1, 360);
SRP_WOA13 = ncread(SRP_path_WOA13, "p_an");
SRP_WOA13_surface = transpose(SRP_WOA13(:, :, 1));
% regrid the SRP 
SRP_WOA13_lon = double(SRP_WOA13_lon);
SRP_WOA13_lat = double(SRP_WOA13_lat);

%% WOA 2023  SRP mapping
% reading data
SRP_path_WOA23 = 'woa23_all_p00_01.nc';
SRP_WOA23_lon = ncread(SRP_path_WOA23,"lon");
SRP_WOA23_lon = repmat(SRP_WOA23_lon', 180, 1);
SRP_WOA23_lat = ncread(SRP_path_WOA23,"lat");
SRP_WOA23_lat = repmat(SRP_WOA23_lat, 1, 360);
SRP_WOA23 = ncread(SRP_path_WOA23, "p_an");
SRP_WOA23_surface = transpose(SRP_WOA23(:, :, 1));
% regrid the SRP 
SRP_WOA23_lon = double(SRP_WOA23_lon);
SRP_WOA23_lat = double(SRP_WOA23_lat);

%% Glodap
SRP_path_Glodap = 'GLODAPv2.2016b.PO4.nc';
SRP_Glodap_lon = ncread(SRP_path_Glodap,"lon");
SRP_Glodap_lon = repmat(SRP_Glodap_lon', 180, 1);
SRP_Glodap_lon(SRP_Glodap_lon > 180) = SRP_Glodap_lon(SRP_Glodap_lon > 180) - 360;
SRP_Glodap_lon = SRP_Glodap_lon(:, [161:360 1:160]);
SRP_Glodap_lat = ncread(SRP_path_Glodap,"lat");
SRP_Glodap_lat = repmat(SRP_Glodap_lat, 1, 360);
SRP_Glodap = ncread(SRP_path_Glodap, "PO4");
SRP_Glodap_surface = transpose(SRP_Glodap(:, :, 1));
SRP_Glodap_surface = SRP_Glodap_surface(:, [161:360 1:160]);
%% in originall Glodap, the Gulf of Mexico, the Caribbean Sea masked.
% replace using WOA2023
SRP_Glodap_surface_plus_WOA23 = SRP_Glodap_surface;
SRP_Glodap_surface_plus_WOA23(isnan(SRP_Glodap_surface)) = SRP_WOA23_surface(isnan(SRP_Glodap_surface));
SRP_Glodap_surface_plus_WOA23(SRP_Glodap_surface_plus_WOA23<0) =SRP_WOA23_surface(SRP_Glodap_surface_plus_WOA23<0);
%%
%writematrix(SRP_WOA13_surface(:), "SRP_WOA13.csv") % for analysis in R
%writematrix(SRP_WOA23_surface(:), "SRP_WOA23.csv") % for analysis in R
%writematrix(SRP_Glodap_surface(:), "SRP_Glodap.csv") % for analysis in R
%writematrix(SRP_Glodap_surface_plus_WOA23(:), "SRP_Glodap_plus_WOA23.csv") % for analysis in R
%%
lat_vector = -89.5:1:89.5;
lon_vector = -179.5:1:179.5;
[lon_mat, lat_mat] = meshgrid(lon_vector, lat_vector);
%% SRP Glodap plus WOA23
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
m_contourf(lon_re, lat_re, SRP_Glodap_surface_plus_WOA23(:, ind),50,'LineStyle','none','EdgeColor','none');
colormap(cmocean('balance'));
%clim([0 0.5]);
%colormap(flipud(hot));
h=colorbar;
set(h, 'FontName', 'Times New Roman', 'FontSize',25);
set(get(h,'ylabel'),'String','DIP (µM)','FontSize',25, 'FontName', 'Times New Roman');
m_coast('patch',[.6 .6 .6]);
m_grid('tickdir','out','linewi',2,'FontSize', 20, 'FontName', 'Times New Roman');
text(-3, 1.6, 'a', 'fontsize',35, 'FontName','Times New Roman');

%%
%print(gcf, 'Figures/fig_DIP', '-dpng', '-r600');
%% DOP mapping
% reading data
load MLpredic.mat % from Liang et al., 2022, Nature geoscience
load OCIM2_CTL_He.mat % from Liang et al., 2022, Nature geoscience
M3d = output.M3d;
grid = output.grid;
DOP_fit = M3d*0+nan;
DOP_fit(1:16380) = mean([boostedtree,SVM,Gaussian],2);
DOP_fit(M3d(:) == 0) = NaN;
DOP_fit_surface_Liang = DOP_fit(:,:,1);
DOP_lat_Liang= grid.YT;
DOP_lon_Liang=grid.XT;
%%
DOP_lon_Liang_rearrange = DOP_lon_Liang(:, [91:180 1:90]);
DOP_lon_Liang_rearrange(DOP_lon_Liang_rearrange>180) = DOP_lon_Liang_rearrange(DOP_lon_Liang_rearrange>180)-360;
DOP_fit_surface_Liang_rearrange = DOP_fit_surface_Liang(:, [91:180 1:90]);
%%
DOP_fit_surface_regrid = interp2(DOP_lon_Liang_rearrange, DOP_lat_Liang, DOP_fit_surface_Liang_rearrange, lon_mat, lat_mat);
% no value at lon +-179.5
%%
lon_missing = repmat([179.5 180.5], 180, 1);
lat_missing = repmat([-89.5:1:89.5]', 1, 2);
DOP_fit_surface_regrid_missing = interp2(DOP_lon_Liang, DOP_lat_Liang, DOP_fit_surface_Liang, lon_missing, lat_missing);
%%
DOP_fit_surface_regrid(:, [1 360]) = DOP_fit_surface_regrid_missing;
%%
% writematrix(DOP_fit_surface_regrid(:), "DOP_surface.csv") % for analysis in R
%% figure DOP check
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
m_contourf(lon_re, lat_re, DOP_fit_surface_regrid(:, ind),50,'LineStyle','none','EdgeColor','none');
colormap(cmocean('balance'));
%clim([0 0.5]);
%colormap(flipud(hot));
h=colorbar;
set(h, 'FontSize',25, 'FontName', 'Times New Roman')
set(get(h,'ylabel'),'String','DOP (µM)','FontSize',25, 'FontName', 'Times New Roman');
m_coast('patch',[.6 .6 .6]);
m_grid('tickdir','out','linewi',2, 'FontSize',20, 'FontName', 'Times New Roman');
text(-3, 1.6, 'b', 'fontsize',35, 'FontName', 'Times New Roman');
%%
%print(gcf, 'Figures/fig_DOP', '-dpng', '-r600');
%% DOP/TDP
DOP_to_TDP_ratio_mat = DOP_fit_surface_regrid./(DOP_fit_surface_regrid + SRP_Glodap_surface_plus_WOA23);
DOP_to_TDP_ratio_array = repmat(DOP_to_TDP_ratio_mat, [1 1 12]);

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
m_contourf(lon_re, lat_re, DOP_to_TDP_ratio_mat(:, ind),50,'LineStyle','none','EdgeColor','none');
colormap(cmocean('balance'));
clim([0 1]);
%colormap(flipud(hot));
h=colorbar;
set(h, 'FontSize',25, 'FontName', 'Times New Roman')
set(get(h,'ylabel'),'String','DOP to TDP ratio','fontsize',25);
m_coast('patch',[.5 .5 .5]);
m_grid('tickdir','out','linewi',2, 'backcolor', [0.95 0.95 0.95], 'FontSize',20, 'FontName', 'Times New Roman');
%text(-3, 1.6, 'd', 'fontsize',25)
text(-3, 1.6, 'c', 'fontsize',35, 'FontName', 'Times New Roman')
%%
%print(gcf, 'Figures/fig_DOP_TDP', '-dpng', '-r600');
