clc;clear;close all;

%% parameters (make sure if "sheet" and "week_num" are correct every time)
% fitting folder
fitting_folder = 'fitting_20210614';
mkdir(fitting_folder);

% fitting wl
wl_res = 1;
target_wl = 450:wl_res:700;
% save(fullfile(fitting_folder, 'fitting_wl.txt'), 'target_wl', '-ascii', '-tabs');

% chromophore file and file path
epsilon_folder = 'epsilon';
epsilon_files = {'epsilon.txt', ...    % HbO2 & Hb in order
                 'bilirubin_extinction.txt', ...
                 'biliverdin_dimethyl_ester_extinction.txt', ...
                 'methemoglobin_extinction.txt'};

% live data file, file path, sheet, subject{sheet} name, and week number,
% and how many wavelengths will shift
liveData_folder = 'skinBruise_liveData';
liveData_file = '20200224_BruiseData.xlsx';
sheet = 3; % only 1 or 3
subject = {'FY','' ,'ZO'}; % correspond to sheet number
week_num = 4; % 5 for FY, 4 for ZO
wl_shift = 0; % 20200807 updated

% each pathlength of wl
pl_file = 'pathlength_0803/pathlength.txt';

%% read and process epsilon data
% read data
epsilon_raw = cell(1, length(epsilon_files)+1);

% read HbO2 & Hb
temp_epsilon = importdata(fullfile(epsilon_folder, epsilon_files{1}));
epsilon_raw{1} = [temp_epsilon(:, 1), temp_epsilon(:, 2)];
epsilon_raw{2} = [temp_epsilon(:, 1), temp_epsilon(:, 3)];
% read bilirubin & biliverdin
for index_epsilon = 3:length(epsilon_files)
    epsilon_raw{index_epsilon} = importdata(fullfile(epsilon_folder, epsilon_files{index_epsilon-1}));
    epsilon_raw{index_epsilon} = epsilon_raw{index_epsilon}.data;
end
% read methemoglobin
epsilon_raw{end} = importdata(fullfile(epsilon_folder, epsilon_files{end}));

% interpolation (get preferred wl epsilon)
target_wl = reshape(target_wl, [], 1);

for index=1:length(epsilon_raw)
    epsilon_interp{index} = [target_wl, interp1(epsilon_raw{index}(:, 1), epsilon_raw{index}(:, 2), target_wl)];
end

% concatenate all epsilon into a matrix
% [HbO2, Hb, bilirubin, biliverdin, methemoglobin]
epsilon_used = [epsilon_interp{1}(:, 2), ...
                epsilon_interp{2}(:, 2), ...
                epsilon_interp{3}(:, 2), ...
                epsilon_interp{4}(:, 2), ...
                epsilon_interp{5}(:, 2)];

%% read and process liveData
% load .xlsx file
liveData_raw.(subject{sheet}) = xlsread(fullfile(liveData_folder, liveData_file), sheet);

% slice data to different weeks, bruise and normal
liveData_categories = {'bruised', 'normal'};
idx = 0;
for index_week = 1:week_num
    for index_category=1:length(liveData_categories)
        liveData_seg.(subject{sheet}).(['week_', num2str(index_week)]).(liveData_categories{index_category}) = ...
                                      [liveData_raw.(subject{sheet})(:, 1 + idx*14), liveData_raw.(subject{sheet})(:, 13 + idx*14)];
        idx = idx+1;
        % do wavelength shift: 20200807 updated
        liveData_seg.(subject{sheet}).(['week_', num2str(index_week)]).(liveData_categories{index_category})(:, 1) = ...
            liveData_seg.(subject{sheet}).(['week_', num2str(index_week)]).(liveData_categories{index_category})(:, 1) - wl_shift;
    end
end

% interpolation (get preferred wl liveData)
for index_week = 1:week_num
    for index_category=1:length(liveData_categories)
        temp_data = liveData_seg.(subject{sheet}).(['week_', num2str(index_week)]).(liveData_categories{index_category});
        liveData_used.(subject{sheet}).(['week_', num2str(index_week)]).(liveData_categories{index_category}) = ...
                                                                   [target_wl, interp1(temp_data(:, 1), temp_data(:, 2), target_wl)];
    end
end

%% read and process pathlength data
% read data
pl_raw = load(pl_file);

% interpolation (get preferred wl pathlength)
% [wl_num, 2]
pl_used = [target_wl, interp1(pl_raw(:, 1), pl_raw(:, 2), target_wl)];

%% main
% calculate delta OD and delta mua
for index_week = 1:week_num
    normal = liveData_used.(subject{sheet}).(['week_', num2str(index_week)]).('normal');
    bruised = liveData_used.(subject{sheet}).(['week_', num2str(index_week)]).('bruised');
    % [wl_num, 2] (no unit)
    liveData_used.(subject{sheet}).(['week_', num2str(index_week)]).('delta_OD') = [target_wl, log(normal(:, 2)./ bruised(:, 2))];
    % [wl_num, 2] (1/cm)
    liveData_used.(subject{sheet}).(['week_', num2str(index_week)]).('delta_mua') = ...
                      [target_wl, liveData_used.(subject{sheet}).(['week_', num2str(index_week)]).('delta_OD')(:, 2)./ pl_used(:, 2)];
end

% calculate fitted delta_c (M) (this study target)
% 2.303 * epsilon * delta_concentration = delta_mua (Ax=b, x=A\b) ------>
% (2.303 * [wl_num, epsilon_species]) * [delta_c, 1] = [delta_mua, 1]
% [delta_c, 1] = (2.303 * [wl_num, epsilon_species]) \ [delta_mua, 1]
for index_week = 1:week_num
    % [chromophore_num, 1] (M)
    liveData_used.(subject{sheet}).(['week_', num2str(index_week)]).('fitted_delta_c') = (2.303 * epsilon_used) \ ...
                                                 liveData_used.(subject{sheet}).(['week_', num2str(index_week)]).('delta_mua')(:, 2);
end

% calculate fitted delta_mua
for index_week = 1:week_num
    % [wl_num, 2] (1/cm)
    liveData_used.(subject{sheet}).(['week_', num2str(index_week)]).('fitted_delta_mua') = ...
           [target_wl, (2.303 * epsilon_used) * liveData_used.(subject{sheet}).(['week_', num2str(index_week)]).('fitted_delta_c')];
%     error = (2.303 * epsilon_used) * liveData_used.(subject{sheet}).(['week_', num2str(index_week)]).('fitted_delta_c') - ...
%                                       liveData_used.(subject{sheet}).(['week_', num2str(index_week)]).('delta_mua')(:, 2);
%     liveData_used.(subject{sheet}).(['week_', num2str(index_week)]).('error') = sum(error.^2);
end

% calculate fitted delta_OD
for index_week = 1:week_num
    % [wl_num, 2] (1/cm)
    liveData_used.(subject{sheet}).(['week_', num2str(index_week)]).('fitted_delta_OD') = ...
           [target_wl, liveData_used.(subject{sheet}).(['week_', num2str(index_week)]).('fitted_delta_mua')(:, 2) .* pl_used(:, 2)];
%     error = (2.303 * epsilon_used) * liveData_used.(subject{sheet}).(['week_', num2str(index_week)]).('fitted_delta_c') - ...
%                                       liveData_used.(subject{sheet}).(['week_', num2str(index_week)]).('delta_mua')(:, 2);
%     liveData_used.(subject{sheet}).(['week_', num2str(index_week)]).('error') = sum(error.^2);
end