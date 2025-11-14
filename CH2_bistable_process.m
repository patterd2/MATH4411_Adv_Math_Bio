function CH2_bistable_process
%% parameters
k1=0.00025;
k2=0.18;
k3=37.5;
k4=1800;

T = 100; % final time

%% Solve ODE and plot solutions
[t1,z1] = ode45(@(t,z) -k1*z(1)*z(1)*z(1)+k2*z(1)*z(1)-k3*z(1)+k4,[0 1],0);
[t2,z2] = ode45(@(t,z) -k1*z(1)*z(1)*z(1)+k2*z(1)*z(1)-k3*z(1)+k4,[0 1],200);
[t3,z3] = ode45(@(t,z) -k1*z(1)*z(1)*z(1)+k2*z(1)*z(1)-k3*z(1)+k4,[0 1],300);
[t4,z4] = ode45(@(t,z) -k1*z(1)*z(1)*z(1)+k2*z(1)*z(1)-k3*z(1)+k4,[0 1],500);

figure(1);
set(gca,'Fontsize',18);
plot(t1,z1,'r','Linewidth',3);
hold on;
plot(t2,z2,'k','Linewidth',3);
plot(t3,z3,'g','Linewidth',3);
plot(t4,z4,'b','Linewidth',3);

% Only uncomment these lines for the fixed points with k4 = 2200
line([0 1],[100 100],'LineStyle',':','Color','k','Linewidth',2);
line([0 1],[220 220],'LineStyle',':','Color','k','Linewidth',2);
line([0 1],[400 400],'LineStyle',':','Color','k','Linewidth',2);

xlabel('time [min]');
ylabel('number of A molecules');
hh=legend('A(0)=0','A(0)=200','A(0)=300','A(0)=500');
set(hh,'location','east','Fontsize',18);
axis([0 1 0 550]);
box on;
set(gca,'Fontsize',20);
grid on;

%% SSA to simulate the stochastic process and plot solutions vs ODE solutions
X=0;
time=0;
timeSSA=0;
kk=0;

while (time<T)
    timeSSAprev=timeSSA;
    timefinSSA=1;
    kk=kk+1;
    while (timeSSA<timeSSAprev+timefinSSA)
        rr=rand(2,1);
        timeSSA=timeSSA+1;
        a0=k1*X*(X-1)*(X-2)+k2*X*(X-1)+k3*X+k4;
        tau=(1/a0)*log(1/rr(1));
        if (rr(2)*a0<(k2*X*(X-1)+k4))
            X=X+1;
        else
            X=X-1;
        end
        time=time+tau;
    end
    XX(kk)=X;
    tt(kk)=time;
end

[t0,z0] = ode45(@(t,z) -k1*z(1)*z(1)*z(1)+k2*z(1)*z(1)-k3*z(1)+k4,[0 T],0);
[t,z] = ode45(@(t,z) -k1*z(1)*z(1)*z(1)+k2*z(1)*z(1)-k3*z(1)+k4,[0 T],500);

figure(2);
line([5 5],[0 0],'Color','b','Linewidth',4);
hold on;
line([5 5],[0 0],'Color','r','Linewidth',4);
h=stairs(tt,XX);
set(h,'Color','b','Linewidth',1);
plot(t,z,'r','Linewidth',3);
plot(t0,z0,'r','Linewidth',3);
xlabel('time [min]');
ylabel('number of A molecules');
hold on;

% Only uncomment these lines for the fixed points with k4 = 2200
line([0 T],[100 100],'LineStyle',':','Color','k','Linewidth',2);
line([0 T],[220 220],'LineStyle',':','Color','k','Linewidth',2);
line([0 T],[400 400],'LineStyle',':','Color','k','Linewidth',2);
 
legend('stochastic','deterministic');
axis([0 T 0 550]);
box on;
set(gca,'Fontsize',20);
grid on;

figure(1);
