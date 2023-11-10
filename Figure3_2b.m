function Figure3_2b

k1=0.001;
k2=0.75;
k3=165;
k4=10000;
k5=200;

ps=zeros(601,1);
xg=0:600;
ps=exp((-3*k1*xg.*xg.*xg.*xg+4*k2*xg.*xg.*xg-6*k3*xg.*xg+12*k4*xg)/(6*k5*k5));
ps=ps/sum(ps);

data_Figure3_2b = load('data_Figure3_2b.dat');

figure(2);
h1=bar([700],[0]);
hold on;
stem([1:600],data_Figure3_2b,'y','Markersize',0.05,'Color',[0.25 0.25 1]);
h2=plot(xg,ps,'r','Linewidth',3);
xlabel('x');
ylabel('p_s(x)');
h3=legend([h1,h2],'SSA','Fokker-Planck');
set(gca,'XTick',[0 100 200 300 400 500]);
axis([0 500 0 0.01]);
set(gca,'Fontsize',20);
grid on;

dx=0.5;
psfun = @(x)exp((-3*k1*x.*x.*x.*x+4*k2*x.*x.*x-6*k3*x.*x+12*k4*x)/(6*k5*k5));
overlinetau=0;
for ii=100.5:0.5:249.5
    xv=ii;
    integral=quad(psfun,10,xv);
    overlinetau=overlinetau+dx*(2/(k5*k5))*(exp(-(-3*k1*xv*xv*xv*xv+4*k2*xv*xv*xv-6*k3*xv*xv+12*k4*xv)/(6*k5*k5)))*integral;
end
    xv=100;
    integral=quad(psfun,10,xv);
    overlinetau=overlinetau+(dx/2)*(2/(k5*k5))*(exp(-(-3*k1*xv*xv*xv*xv+4*k2*xv*xv*xv-6*k3*xv*xv+12*k4*xv)/(6*k5*k5)))*integral;
    xv=250;
    integral=quad(psfun,10,xv);
    overlinetau=overlinetau+(dx/2)*(2/(k5*k5))*(exp(-(-3*k1*xv*xv*xv*xv+4*k2*xv*xv*xv-6*k3*xv*xv+12*k4*xv)/(6*k5*k5)))*integral;

disp(['Switching Time:',num2str(overlinetau)]);

