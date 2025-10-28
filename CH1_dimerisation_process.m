function CH1_dimerisation_process

k1 = 0.0075;
k2 = 1;
T = 100; % max time for each simulation

%% Direct simulations of the dimerisation process
Ainitial = 15;
numberofrealisations = 1000;
p=zeros(24,1);

for i=1:numberofrealisations
    A=Ainitial;
    time=0;
    k=1;
    Aplot(1,i)=A;
    timeplot(1,i)=0;
    while (time < T)
        rr=rand(2,1);
        a0=k1*A*(A-1)+k2;
        time = time+(1/a0)*log(1/rr(1));
        if (rr(2)*a0<k1*A*(A-1))
            A=A-2;
        else
            A=A+1;
        end
        k=k+1;
        Aplot(k,i)=A;
        timeplot(k,i)=time;
    end
    if (A<23.5)
        p(A+1)=p(A+1)+1;
    end
end
p=p/numberofrealisations;
% p(n+1) is the ESTIMATED stationary probability that A=n for n=0,1,2,3,.....

%% Use the formula for the PGF to calculate the stationary mean and variance
x= -1:0.01:1;
w=size(x);
for i=1:w(2)
    Gs(i)=sqrt(1+x(i))*besseli(1,2*sqrt((1+x(i))*k2/k1));
    Gsder(i)=1/(2*sqrt(1+x(i)))*besseli(1,2*sqrt((1+x(i))*k2/k1))+sqrt(k2/k1)/2*(besseli(2,2*sqrt((1+x(i))*k2/k1))+besseli(0,2*sqrt((1+x(i))*k2/k1)));
    Gsderder(i)=(k2/k1)*(1/(1+x(i)))*Gs(i);
end
Gs=Gs/(sqrt(2)*besseli(1,2*sqrt(2*k2/k1)));
Gsder=Gsder/(sqrt(2)*besseli(1,2*sqrt(2*k2/k1)));
Gsderder=Gsderder/(sqrt(2)*besseli(1,2*sqrt(2*k2/k1)));

% Calculate the mean and variance from the PGF
% i.e. evaluate the appropriate derivatives at x = 1
Ms = Gsder(w(2));
Vs = Gsderder(w(2))+Ms-Ms*Ms;
%% Calculate the stationary distribution using the analytic formula
% This formula is derived in Problem Class 2, Question 3
for n = 0:25
    phi(n+1) = ((sqrt(k2/k1)^n)/factorial(n))*besseli(n-1,2*(sqrt(k2/k1)));
end
phi=phi/sum(phi);

%% Use a built in ODE solver to solve the law of mass action ODE
[tdet,Adet] = ode45(@(t,z) -2*k1*z*z+k2,[0 T], Ainitial);

%% Plotting of sample paths + solution to ODE + stationary mean & variance
figure;
plot(tdet,Adet,'--k','Linewidth',4);
hold on
plot(tdet,Ms*ones(1,length(tdet)),'--r','Linewidth',5);
plot(tdet,(Ms+2*sqrt(Vs))*ones(1,length(tdet)),'-b','Linewidth',5);
plot(tdet,(Ms-2*sqrt(Vs))*ones(1,length(tdet)),'-b','Linewidth',5);
for i = 1:min(25,numberofrealisations)
    [~, t2] = max(timeplot(:,i));
    h=stairs(timeplot(1:t2,i),Aplot(1:t2,i));
    set(h,'Linewidth',1);
    hold on;
end
plot(tdet,Adet,'--k','Linewidth',5);
hold on;
plot(tdet,Ms*ones(1,length(tdet)),'--r','Linewidth',5);
plot(tdet,(Ms+2*sqrt(Vs))*ones(1,length(tdet)),'-b','Linewidth',4);
plot(tdet,(Ms-2*sqrt(Vs))*ones(1,length(tdet)),'-b','Linewidth',4);
xlabel('time [sec]');
ylabel('number of molecules');
legend('mass action model','M_s','M_s \pm 2 (V_s)^{1/2}','location','southwest');
axis tight;
xlim([0 T]);
ylim([0 max(max(Aplot))+1]);
set(gca,'Fontsize',20);
grid on;
%keyboard;
%% Plotting of the empirical and stationary PMFs
figure;
h1=bar(0:23,p);
hold on;
h2=plot(0:25,phi,'r','Linewidth',4);
xlabel('number of molecules');
ylabel('stationary distribution');
hh=legend([h1 h2],'Gillespie SSA','analysis');
axis tight;
xlim([0 25]);
set(gca,'Fontsize',20);
grid on;
