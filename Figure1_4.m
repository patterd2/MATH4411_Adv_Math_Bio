function Figure1_4

% Matlab code for Figure 1.4 in book:
% Radek Erban and Jon Chapman, Stochastic Modelling of Reaction-Diffusion
% Processes, Cambridge University Press, 2019

close all;

z=[0:0.02:4];

figure(1);
set(gca,'Fontsize',18);
plot(z,besseli(1,z),'r','Linewidth',3);
hold on;
plot(z,besselk(1,z),'k-.','Linewidth',3);
xlabel('$z$','interpreter','latex');
ylabel('modified Bessel functions','interpreter','latex');
hh=legend('$I_1(z)$','$K_1(z)$');
set(hh,'interpreter','latex','location','north','Fontsize',18);
set(gca,'XTick',[0 0.5 1 1.5 2 2.5 3]);
axis([0 3.3 0 5]);
set(gca,'Fontsize',18);

k1=0.005;
k2=1;
x=[0:0.01:1];
w=size(x);
for i=1:w(2)
    Gs(i)=sqrt(1+x(i))*besseli(1,2*sqrt((1+x(i))*k2/k1));
    Gsder(i)=1/(2*sqrt(1+x(i)))*besseli(1,2*sqrt((1+x(i))*k2/k1))+sqrt(k2/k1)/2*(besseli(2,2*sqrt((1+x(i))*k2/k1))+besseli(0,2*sqrt((1+x(i))*k2/k1)));
    Gsderder(i)=(k2/k1)*(1/(1+x(i)))*Gs(i);
end;
Gs=Gs/(sqrt(2)*besseli(1,2*sqrt(2*k2/k1)));
Gsder=Gsder/(sqrt(2)*besseli(1,2*sqrt(2*k2/k1)));
Gsderder=Gsderder/(sqrt(2)*besseli(1,2*sqrt(2*k2/k1)));

Ms=Gsder(w(2));
disp(['M_s = ',num2str(Ms)]);
Vs=Gsderder(w(2))+Ms-Ms*Ms;
disp(['V_s = ',num2str(Vs)]);
disp(['square root of Vs = ',num2str(sqrt(Vs))]);

figure(2);
set(gca,'Fontsize',18);
h1=line([10 10],[10 10],'Color','r');
hold on;
h2=line([10 10],[10 10],'Color','b');
[AX,h3,h4] = plotyy(x,Gs,x,Gsder);
xlabel('$x$','interpreter','latex');
set(AX(1),'Xlim',[0.5 1]); 
set(AX(2),'Xlim',[0.5 1]); 
set(AX(1),'Ylim',[0 1]); 
set(AX(2),'Ylim',[0 15]); 
set(AX(1),'Fontsize',20);
set(AX(2),'Fontsize',20);
set(get(AX(1),'Ylabel'),'String','$G_s(x)$','Fontsize',18,'interpreter','latex'); 
set(get(AX(2),'Ylabel'),'String','${\rm d}G_s/{\rm d}x (x)$','Fontsize',18,'interpreter','latex'); 
set(h1,'Linewidth',3);
set(h2,'Linewidth',3);
set(h2,'Linestyle','--');
set(h3,'Linewidth',3);
set(h3,'Color','r');
set(h4,'Linewidth',3);
set(h4,'Linestyle','--');
set(h4,'Color','b');
set(AX(1),'YColor','r');
set(AX(2),'YColor','b');
set(get(AX(2),'Xlabel'),'String','$x$','Rotation',0,'Fontsize',18,'interpreter','latex'); 
xlabel('$x$','interpreter','latex');
set(AX(1),'XTick',[0.5 0.6 0.7 0.8 0.9 1]);
set(AX(1),'YTick',[0 0.2 0.4 0.6 0.8 1]);
set(AX(2),'YTick',[0 1 2 3 4 5 6 7 8 9 10 11 12]);
set(AX(2),'YTicklabel',{'0','','2','','4','','6','','8','','10','',''});
hh=legend('$G_s$','${\rm d}G_s/{\rm d}x$');
set(hh,'interpreter','latex','location','northwest','Fontsize',18);
set(gca,'Fontsize',18);


