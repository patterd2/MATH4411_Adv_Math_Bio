function CH4_compartment_RD

D=0.0001;
h=0.01;
xdet=[h/2:h:1-h/2];
w=size(xdet);
boxdeterm=w(2);
dt=0.01;
Diff=D/(h*h);
detden=zeros(boxdeterm,1);

k2=0.3125;
k1=0.001;

for t=0:dt:10*60
    detdenold=detden;
    ii=1;
    detden(ii)=detdenold(ii)+dt*Diff*(detdenold(ii+1)-detdenold(ii))-dt*k1*detdenold(ii);
    for ii=2:boxdeterm-1
        detden(ii)=detdenold(ii)+dt*Diff*(detdenold(ii+1)+detdenold(ii-1)-2*detdenold(ii))-dt*k1*detdenold(ii);
    end
    ii=boxdeterm;
    detden(ii)=detdenold(ii)+dt*Diff*(detdenold(ii-1)-detdenold(ii))-dt*k1*detdenold(ii);

    for ii=1:boxdeterm/5
        detden(ii)=detden(ii)+dt*k2;
    end
    % figure(1);
    % bar(mesh,histogram(:,1));
    % hold on;
    % set(gca,'Fontsize',20);
    % grid on;
end

histogram=load('data_Figure6_1.dat');
ww=size(histogram);
boxes=ww(1);
mesh=[1/(2*boxes):1/boxes:1-1/(2*boxes)]';

figure;
plot([0 xdet 1],[detden(1); detden; detden(boxdeterm)],'r','Linewidth',4);
hold on;
hb1=bar(mesh,histogram(:,1),'b');
plot([0 xdet 1],[detden(1); detden; detden(boxdeterm)],'r','Linewidth',4);
box on;
xlabel('x','interpreter','tex');
ylabel('number of molecules','interpreter','tex');
text(0.15,142,'time=10 min','Fontsize',20,'interpreter','tex');
set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1]);
hh1=legend('PDE approx.','Compartmental model');
axis([0 1 0 155]);
set(gca,'Fontsize',20);
grid on;

for t=10*60+dt:dt:30*60
    detdenold=detden;
    ii=1;
    detden(ii)=detdenold(ii)+dt*Diff*(detdenold(ii+1)-detdenold(ii))-dt*k1*detdenold(ii);
    for ii=2:boxdeterm-1
        detden(ii)=detdenold(ii)+dt*Diff*(detdenold(ii+1)+detdenold(ii-1)-2*detdenold(ii))-dt*k1*detdenold(ii);
    end
    ii=boxdeterm;
    detden(ii)=detdenold(ii)+dt*Diff*(detdenold(ii-1)-detdenold(ii))-dt*k1*detdenold(ii);
    for ii=1:boxdeterm/5
        detden(ii)=detden(ii)+dt*k2;
    end
end

figure;
plot([0 xdet 1],[detden(1); detden; detden(boxdeterm)],'r','Linewidth',4);
hold on;
hb2=bar(mesh,histogram(:,2),'b');
plot([0 xdet 1],[detden(1); detden; detden(boxdeterm)],'r','Linewidth',4);
box on;
xlabel('x','interpreter','tex');
ylabel('number of molecules','interpreter','tex');
text(0.15,142,'time=30 min','Fontsize',20,'interpreter','tex');
hh2=legend('PDE approx.','Compartmental model');
set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1]);
set(gca,'YTick',[0 20 40 60 80 100 120 140]);
axis([0 1 0 155]);
set(gca,'Fontsize',20);
grid on;

