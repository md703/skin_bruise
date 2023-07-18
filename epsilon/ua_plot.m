%
w = 250:2:1000;
ua_water = load('water_mua.txt');       ua_water = ua_water(:,2);
ua_mel = load('mel_mua.txt');           ua_mel = ua_mel(:,2);
a = load('epsilon.txt');                molar_hc = a(:,2); molar_dehc = a(:,3);
ua_up = load('epsilon_2.txt');          ua_up = ua_up(:,2);
ua_colla = load('collagen_mua_data.txt');ua_colla = ua_colla(:,2);

%%%
hc = 0.4;
sto2 = 0.8;
colla = 0.2;
%¤U¼h§l¦¬¦±½u without water
ua_bt = hc*(sto2*2.303*molar_hc/64532+(1-sto2)*2.303*molar_dehc/64500) + colla*ua_colla;
% with water
ua_bt_water = ua_bt + (1- (hc/150 + colla))*ua_water;
figure();
plot(w,ua_bt,'b',w,ua_bt_water,'r');

figure();
percent = (ua_bt_water-ua_bt)./ua_bt;
plot(w,percent)
ylabel('100%')
title('percentage add')