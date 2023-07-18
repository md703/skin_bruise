eu = load('eumelanin.txt');
ph = load('pheomelanin.txt');
w = 400:2:800;
w_eu = eu(:,1);
w_ph = ph(:,1);
eu_ = interp1(w_eu,eu(:,2),w);
ph_ = interp1(w_ph,ph(:,2),w);
skin = (1.70e12)*(w.^(-3.48));
%{
eu_n = eu_./eu_(101);
ph_n = ph_./ph_(101);
skin_n = skin./skin(101);
figure();
plot(w,eu_n,'b',w,ph_n,'g',w,skin_n,'r')
%}

xdata = w;
ydata = skin;
x0 = [9 1 100];
fun = @(x,xdata)x(3).*((x(1)/(sum(x))).*(eu_)+(x(2)./sum((x))).*(ph_));
x = lsqcurvefit(fun,x0,xdata,ydata)

sim = x(3).*((x(1)/(sum(x))).*(eu_)+(x(2)./sum((x))).*(ph_));
figure();
plot(xdata,ydata,'r',xdata,sim,'m')
hold on
%plot(w,eu_n,'b',w,ph_n,'g')
hold off


