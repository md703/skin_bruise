%in vivo, skin, melanin
%from: omlc.rog/spectra

w = 250:2:1000
mua_mel = (1.7e12)*(w.^-3.48);

figure();
plot(w,mua_mel)
figure();
semilogy(w,mua_mel)

list = [w' mua_mel']