function CH1_Gillespie_example

k1=0.001;
k2=0.01;
k3=1.2;
k4=1;

Ainitial=15;
Binitial=15;
numberofrealisations=5;
T = 50;

for i = 1:numberofrealisations
    A=Ainitial;
    B=Binitial;
    time=0;
    k=1;
    Aplot(1,i)=A; 
    Bplot(1,i)=B;
    timeplot(1,i)=0;
    while (time<T)
        rr=rand(2,1);
        a0=k1*A*(A-1)+k2*A*B+k3+k4;
        time=time+(1/a0)*log(1/rr(1));
        ss=k1*A*(A-1);
        if (rr(2)*a0<ss)
            A=A-2;
        else
            ss=ss+k2*A*B;
            if (rr(2)*a0<ss)
                A=A-1;
                B=B-1;
            else
                ss=ss+k3;
                if (rr(2)*a0<ss)
                    A=A+1;
                else
                    B=B+1;
                end
            end
        end
        k=k+1;
        Aplot(k,i)=A;
        Bplot(k,i)=B;
        timeplot(k,i)=time;
    end
end

[tdet,det] = ode45(@myode,[0 T],[Ainitial;Binitial]);

%% Plotting
figure;
set(gca,'Fontsize',18);
plot(tdet,det(:,1),'--k','Linewidth',4);
hold on;
for j = 1:numberofrealisations
    [~, t2] = max(timeplot(:,j));
    stairs(timeplot(1:t2,j),Aplot(1:t2,j),'Linewidth',2);
end
plot(tdet,det(:,1),'--k','Linewidth',4);
xlabel('time [sec]');
ylabel('number of A molecules');
legend('solution of ODEs');
axis tight;
ylim([0 25]);
set(gca,'Fontsize',20);
grid on;

figure;
plot(tdet,det(:,2),'--k','Linewidth',4);
hold on;
for j = 1:numberofrealisations
    [~, t2] = max(timeplot(:,j));
    stairs(timeplot(1:t2,j),Bplot(1:t2,j),'Linewidth',2);
end
plot(tdet,det(:,2),'--k','Linewidth',4);
xlabel('time [sec]');
ylabel('number of B molecules');
legend('solution of ODEs');
axis tight
ylim([0 25]);
set(gca,'Fontsize',20);
grid on;

% figure;
% for j = 1:numberofrealisations
%     [~, t2] = max(timeplot(:,j));
%     scatter(Aplot(1:t2,j),Bplot(1:t2,j));
% end
% xlabel('number of A molecules');
% ylabel('number of B molecules');
% axis tight
% set(gca,'Fontsize',20);
% grid on;

%% Define RHS of the ODE
function dydt = myode(t,z)
k1=0.001;
k2=0.01;
k3=1.2;
k4=1;
dydt = [k3-2*k1*z(1)*z(1)-k2*z(1)*z(2); k4-k2*z(1)*z(2)];

