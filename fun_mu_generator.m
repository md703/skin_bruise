%{
Genetate the optical parameters from the given tissue parameters

n: number of wavelength
L: number of layers

Inputs:
para: the setting of the tissue parameters
lambda: the wavelength to calculate, a [n*1] array

Outputs:
mu: the optical parameters of each layer, a [n * 2L] array

Benjamin Kao
Last update: 2020/03/03
%}

function mu = fun_mu_generator(para,lambda)

    %% param
    input_dir='epsilon';
    
    vol_water = 0.7;
    
    % mua of         upper layer     Hb & HbO      collagen                water           melanin
    OP_filename_arr={'epsilon_2.txt','epsilon.txt','collagen_mua_data.txt','water_mua.txt','mel_mua.txt'};

    %% load files
    OP_spec_arr=cell(1,length(OP_filename_arr));
    for i=1:length(OP_filename_arr)
        OP_spec_arr{i}=load(fullfile(input_dir,OP_filename_arr{i}));
    end

    %% ckeck the wavelength range
    OPs_wl=[0,inf]; % the min and the max range of the absorber optical parameters
    for i=1:length(OP_filename_arr)
        if OPs_wl(1)<min(OP_spec_arr{i}(:,1))
            OPs_wl(1)=min(OP_spec_arr{i}(:,1));
        end
        if OPs_wl(2)>max(OP_spec_arr{i}(:,1))
            OPs_wl(2)=max(OP_spec_arr{i}(:,1));
        end
    end
    
    assert(min(lambda)>=OPs_wl(1),'ERROR: lambda < literature optiocal parameters');
    assert(max(lambda)<=OPs_wl(2),'ERROR: lambda > literature optiocal parameters');
    
    
    %% main
    
    %% interp
    lambda=reshape(lambda,[],1);
    for i=1:length(OP_filename_arr)
        OP_spec_arr{i}=interp1(OP_spec_arr{i}(:,1),OP_spec_arr{i}(:,2:end),lambda);
    end
    
    %% calculate the OPs
    mu=zeros(length(lambda),para.num_layer*2);
    for L=1:para.num_layer
        %% mua
        if L==1
            mu(:,1)=OP_spec_arr{1};
        else
            ua_hc=((para.tissue_param{L}.sto2*2.303.*OP_spec_arr{2}(:,1)/64535)+((1-para.tissue_param{L}.sto2)*2.303*OP_spec_arr{2}(:,2)/64500))*150;
            vol_colla = 1-vol_water-para.tissue_param{L}.hc-para.tissue_param{L}.mel;
            mu(:,2*L-1)=ua_hc*para.tissue_param{L}.hc+OP_spec_arr{3}*vol_colla+OP_spec_arr{4}*vol_water+OP_spec_arr{5}*para.tissue_param{L}.mel;
        end
        
        %% mus
        if L==1
            mu(:,2)=(1.3*para.tissue_param{2}.A*1000*lambda.^-para.tissue_param{2}.K)/(1-para.tissue_param{1}.g);
        else
            mu(:,2*L)=(para.tissue_param{L}.A*1000*lambda.^-para.tissue_param{L}.K)/(1-para.tissue_param{L}.g);
        end
    end
end