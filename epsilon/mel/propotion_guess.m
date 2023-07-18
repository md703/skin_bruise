%fit melanin eu:ph ¤ñ¨Ò
w = 400:2:800;
eu = load('eumelanin.txt');
ph = load('pheomelanin.txt');
w_eu = eu(:,1);
w_ph = ph(:,1);

eu_ = interp1(w_eu,eu(:,2),w);
ph_ = interp1(w_ph,ph(:,2),w);
skin = (1.70e12)*(w.^(-3.48));
propo = [9 1];
sim_skin = eu_.*(propo(1)/sum(propo)) + ph_.*(propo(2)/sum(propo));

figure();
hold on
%600nm normalize
plot(w,eu_./eu_(101),'b',w,ph_./ph_(101),'g',w,skin./skin(101),'r')
plot(w,sim_skin./sim_skin(101),'--r')
hold off
