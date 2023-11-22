function Figure3_5

k1=0.00025;
k2=0.18;
k3=37.5;
k4=2200;

F = @(x)(k2*x.*(x-1)+k4)./(k1*x.*(x-1).*(x-2)+k2*x.*(x-1)+k3*x+k4);

for ii=1:600
    a(ii)=ii;
    diff2(ii)=k1*a(ii).*(a(ii)-1).*(a(ii)-2)+k2*a(ii).*(a(ii)-1)+k3*a(ii)+k4;
    statdistr(ii)=exp(-(2*a(ii)-4*quad(F,0,a(ii))+log(diff2(ii))+44.4));
end

statdistr=statdistr/sum(statdistr);
potential=-log(statdistr);

figure(1);
set(gca,'Fontsize',18);
h1=bar([700],[0]);
set(h1,'FaceColor','y');
hold on;
h2=plot(a,statdistr,'b','Linewidth',3);

load data_Figure3_5a.dat;
schloglhistogram=data_Figure3_5a;

stem([1:600],schloglhistogram,'y','Markersize',0.01);
plot(a,statdistr,'b','Linewidth',3);
axis([0 600 0 0.015]);
set(gca,'XTick',[0 100 200 300 400 500 600]);
set(gca,'YTick',[0 0.002 0.004 0.006 0.008 0.01 0.012 0.014]);
set(gca,'Fontsize',18);

xlabel('$x$','interpreter','latex');
ylabel('stationary distribution','interpreter','latex');
h3=legend([h1 h2],'$\;$stochastic simulation','$\;$chemical Fokker-Planck');
set(h3,'interpreter','latex','Fontsize',18,'location','northeast');

load data_Figure3_5b.dat;
initialx=data_Figure3_5b(:,1);
exittimes=data_Figure3_5b(:,2);

intp(1)=0;
for ii=2:235
    intp(ii)=intp(ii-1)+statdistr(ii-1)/2+statdistr(ii)/2;
end
tau(235)=0;
for ii=234:-1:1
    tau(ii)=tau(ii+1)+intp(ii+1)/statdistr(ii+1)/(diff2(ii+1))+intp(ii)/statdistr(ii)/(diff2(ii));
end

figure(2);
set(gca,'Fontsize',18);
h4=bar([700],[0]);
set(h4,'FaceColor','y');
hold on;
axis([0 238 0 17]);
set(gca,'XTick',[0 50 100 150 200]);
set(gca,'YTick',[0 2 4 6 8 10 12 14 16]);
xlabel('$y$','interpreter','latex');
ylabel('$\tau(y)$ [min]','interpreter','latex');
h5=plot([0:235],[tau(1) tau],'b','Linewidth',3);
stem(initialx,exittimes,'y','Markersize',0.01);
plot([0:235],[tau(1) tau],'b','Linewidth',3);
h6=legend([h4 h5],'$\;$stochastic simulation','$\;$integral formula');
set(h6,'interpreter','latex','Fontsize',18,'location','southwest');
set(gca,'Fontsize',18);

zz=a(235);
d2xu=diff2(235);
derd2xu=k1*((zz-1)*(zz-2)+zz*(zz-1)+zz*(zz-2))+k2*zz+k2*(zz-1)+k3;
derderphixu=((-4*k2*zz-4*k2*(zz-1)+2*k2+k1*(6*zz-6))*d2xu-(-4*k2*zz*(zz-1)-4*k4+derd2xu)*derd2xu)/(d2xu*d2xu);
zz=a(95);
d2xf=diff2(95);
derd2xf=k1*((zz-1)*(zz-2)+zz*(zz-1)+zz*(zz-2))+k2*zz+k2*(zz-1)+k3;
derderphixf=((-4*k2*zz-4*k2*(zz-1)+2*k2+k1*(6*zz-6))*d2xf-(-4*k2*zz*(zz-1)-4*k4+derd2xf)*derd2xf)/(d2xf*d2xf);
Kramersexittime=3.14159265358979*exp(potential(235)-potential(95))/sqrt(derderphixf*(-derderphixu));
Kramersexittime=2*Kramersexittime/(d2xu);

disp('Mean Exit Time');
disp(['Exact Formula:',num2str(tau(95))]);
disp(['Kramers Approximation:',num2str(Kramersexittime)]);
