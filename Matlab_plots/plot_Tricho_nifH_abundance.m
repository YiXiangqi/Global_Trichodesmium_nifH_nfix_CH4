clear;
clc;
%% set working directory

%% read nifH abundance yearly average data
% `NA` in csv file exported from R should be first replaced with `NaN` 
mlp_Tricho_nifH_abundance_prediction_yearly_average_summary = ...% from Prediction_and_Uncertainty_analysis.Rmd 
    readtable("mlp_Tricho_nifH_abundance_prediction_yearly_average_summary.csv"); 

Tricho_nifH_abundance_pred_mean_mat = ...
    reshape(mlp_Tricho_nifH_abundance_prediction_yearly_average_summary.mean, [180, 360]);
Tricho_nifH_abundance_pred_sd_mat = ...
    reshape(mlp_Tricho_nifH_abundance_prediction_yearly_average_summary.sd, [180, 360]);
%% read Tricho_nifH_abundance_data.csv
Tricho_nifH_abundance_observed_data = readtable("Tricho_nifH_abundance_data.csv");
Tricho_nifH_abundance_observed_data_train = ...
    Tricho_nifH_abundance_observed_data(Tricho_nifH_abundance_observed_data.data_type == "train", :);
Tricho_nifH_abundance_observed_data_test = ...
    Tricho_nifH_abundance_observed_data(Tricho_nifH_abundance_observed_data.data_type == "test", :);


%%
lat_vector = -89.5:1:89.5;
lon_vector = -179.5:1:179.5;
[lon_mat, lat_mat] = meshgrid(lon_vector, lat_vector);

%% figure nifH abundance yearly average:mean 
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

% Adjust observed lon data
lon_nifH_oberved_train = Tricho_nifH_abundance_observed_data_train.lon_regrid;
lon_nifH_oberved_test = Tricho_nifH_abundance_observed_data_test.lon_regrid;
lon_nifH_oberved_train(lon_nifH_oberved_train<0) = lon_nifH_oberved_train(lon_nifH_oberved_train<0)+360;
lon_nifH_oberved_train(lon_nifH_oberved_train>30) = lon_nifH_oberved_train(lon_nifH_oberved_train>30)-360;
lon_nifH_oberved_test(lon_nifH_oberved_test<0) = lon_nifH_oberved_test(lon_nifH_oberved_test<0)+360;
lon_nifH_oberved_test(lon_nifH_oberved_test>30) = lon_nifH_oberved_test(lon_nifH_oberved_test>30)-360;

hold on
m_proj('robinson','lon',[-330 30]); % Specify map domain

% Contour plot
m_contourf(lon_re, lat_re, Tricho_nifH_abundance_pred_mean_mat(:, ind), 50, ...
    'LineStyle', 'none', 'EdgeColor', 'none');

% Scatter plot - TEST data (circle markers)
h1 = m_scatter(lon_nifH_oberved_test, Tricho_nifH_abundance_observed_data_test.lat_regrid, ...
    [], Tricho_nifH_abundance_observed_data_test.Log_Trichodesmium_nifH_integral_regrid, ...
    'filled', 'o', 'MarkerEdgeColor', 'k');

% Scatter plot - TRAINING data (triangle markers)
h2 = m_scatter(lon_nifH_oberved_train, Tricho_nifH_abundance_observed_data_train.lat_regrid, ...
    [], Tricho_nifH_abundance_observed_data_train.Log_Trichodesmium_nifH_integral_regrid, ...
    'filled', '^', 'MarkerEdgeColor', 'k');

% Colormap and colorbar
colormap(flipud(pink(8)))
clim([4 12]);
h = colorbar;
set(get(h, 'ylabel'), ...
    'String', ["{\it Trichodesmium nifH} genen abundance", "(log_{10} copies m^{-2})"], ...
    'fontsize', 14,'FontName','Times New Roman');
set(h, 'FontName','Times New Roman');
% Coastline and grid
m_coast('patch', [.5 .5 .5]);
m_grid('box', 'fancy', 'linestyle', ':', 'gridcolor', 'black', 'backcolor', [0.9 0.9 0.9], ...
    'FontSize', 12, 'FontName','Times New Roman');

% Add legend with horizontal layout
legend_handle = legend([h1, h2], {'Testing Data', 'Training Data'}, ...
    'FontSize', 14, 'Orientation', 'horizontal','FontName','Times New Roman'); % Set orientation to horizontal

% Adjust the legend position manually (bottom-left corner with better alignment)
legend_handle.Position = [0.33, 0.1, 0.2, 0.05]; % [left, bottom, width, height]
text(-3.3, 1.6, 'a', 'FontSize', 20, 'FontWeight','bold', 'FontName','Times New Roman');

%%
%print(gcf, 'Figures/fig_3a', '-dpng', '-r600');
%% figure nifH abundance yearly average:sd 
%tightfig;
% Set figure size
figure
figure_handle = gcf; % Get current figure handle
set(figure_handle, 'Position', [100, 100, 900, 500]); % [left, bottom, width, height]

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

% Contour plot
m_contourf(lon_re, lat_re, Tricho_nifH_abundance_pred_sd_mat(:, ind), 50, ...
    'LineStyle', 'none', 'EdgeColor', 'none');

% Colormap and colorbar
colormap(flipud(pink(9)))
clim([0 1.8]);
h = colorbar;
set(get(h, 'ylabel'), 'String', ...
    ["SD of predicted {\it nifH} abundance", "(log_{10} copies m^{-2})"], 'fontsize', 20);
set(h,'ytick',0:0.2:1.8,'yticklabel',0:0.2:1.8,...
        'tickdir','out', 'FontSize',25, 'FontName', 'Times New Roman');
% Coastline and grid
m_coast('patch', [.5 .5 .5]);
m_grid('box', 'fancy', 'linestyle', ':', 'gridcolor', 'black', 'backcolor', [0.9 0.9 0.9], ...
    'FontSize',20, 'FontName', 'Times New Roman');

%%
%print(gcf, 'Figures/fig_nifH_SD', '-dpng', '-r600');