%{
Simulate the skin reflectance spectrum

Benjamin Kao
Last update: 2020/03/04
%}

clc;clear;close all;

%% param
num_layer=3;
sim_wl=400:100:600;
output_dir='simulation';

n=1.4;

layer_thick=[0.06 0.06]; % in cm, the thickness of layer 1 ~ (L-1), the bottom layer are semi-infinite

% layer 1
% use literature stratum corneum mua
% use musp=1.3 * layer 2 musp
param.tissue_param{1}.g=0.835; % anisotropic factor

% layer 2
param.tissue_param{2}.hc=0; % hemoglobin volume fraction
param.tissue_param{2}.sto2=0; % oxygen saturation
param.tissue_param{2}.mel=0.09; % melanin volume fraction
param.tissue_param{2}.A=1980000; % musp = 1000 * A * lambda ^ -K
param.tissue_param{2}.K=2.833;
param.tissue_param{2}.g=0.75; % anisotropic factor
% layer 3
param.tissue_param{3}.hc=0.05;
param.tissue_param{3}.sto2=0.7;
param.tissue_param{3}.mel=0;
param.tissue_param{3}.A=350322;
param.tissue_param{3}.K=2.631;
param.tissue_param{3}.g=0.715;

%% init
param.num_layer=num_layer;

mkdir(output_dir);
save(fullfile(output_dir,'sim_wl.txt'),'sim_wl','-ascii','-tabs');

%% get the optical parameters
mu=fun_mu_generator(param,sim_wl);

%% start simulation
main_dir=pwd;
for wl=1:length(sim_wl)
    mkdir(fullfile(main_dir,output_dir,['wl_' num2str(wl)]));
    cd(fullfile(main_dir,output_dir,['wl_' num2str(wl)]));
    
    % make input of MCML 
    temp_input=[];
    for L=1:num_layer
        if L<num_layer
            temp_input=[temp_input layer_thick(L) mu(wl,(L*2-1):(L*2)) n param.tissue_param{L}.g];
        else
            temp_input=[temp_input mu(wl,(L*2-1):(L*2)) n param.tissue_param{L}.g];
        end
    end
    save('GPUMC_input.txt','temp_input','-ascii','-tabs');
    
    % send file into MCML_GPU and run simulation
    if isunix
        [~,~]=system('../../bin/MCML_GPU ../../sim_setup.json GPUMC_input.txt GPUMC_output.txt -R -AP');
    elseif ispc
        [~,~]=system('..\..\bin\MCML_GPU.exe ..\..\sim_setup.json GPUMC_input.txt GPUMC_output.txt -R -AP');
    end 
end