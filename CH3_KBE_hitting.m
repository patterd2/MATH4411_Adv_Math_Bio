function CH3_KBE_hitting

k1=0.001;
k2=0.75;
k3=165;
k4=10000;
k5=200;

dt=0.0001;
X(1)=245;
time(1)=0;
steps = 1000;

for i=1:steps
    X(i+1) = X(i)+dt*(-k1*X(i)*X(i)*X(i)+k2*X(i)*X(i)-k3*X(i)+k4)+k5*sqrt(dt)*randn(1,1);
    time(i+1)=i*dt;
end

Xbig=X(1:100:steps+1);
timebig=time(1:100:steps+1);

%% Plotting
figure(1);
plot(time,X,'b','Linewidth',2);
hold on;
h1=plot([1],[1],'b','Linewidth',3);
h2=plot(timebig,Xbig,'r','Linewidth',3);
line([0 steps*dt],[250 250],'Color','k','Linestyle','--','Linewidth',3);
hh=legend([h1 h2],'\Delta_t=10^{-5}','\Delta_t=10^{-3}');
xlim([0 steps*dt]);
ylim([min(min(X),min(Xbig)) max(max(X),max(Xbig))]);
xlabel('t');
ylabel('X(t)');
set(get(gca,'ylabel'),'rotation',0);
set(gca,'Fontsize',20);
grid on;


load data_Figure3_3b_1.dat;
load data_Figure3_3b_2.dat;

initialpoint=data_Figure3_3b_1(:,1);
exittime=data_Figure3_3b_1(:,2);
y=data_Figure3_3b_2(:,1);
tau=data_Figure3_3b_2(:,2);

figure(2);
bar([700],[0]);
hold on;
axis([0 250 0 70]);
set(gca,'XTick',[0 50 100 150 200 250]);
set(gca,'YTick',[0 10 20 30 40 50 60 70]);
xlabel('y');
ylabel('E[\tau (y)]');
set(get(gca,'ylabel'),'rotation',0);
plot(y,tau,'r','Linewidth',3);
stem(initialpoint,exittime,'b','Markersize',0.01);
plot(y,tau,'r','Linewidth',3);
hh=legend('stochastic simulation','integral formula');
set(gca,'Fontsize',20);

