function PS1_Q1_SIS_model

N = 40; % population size
alpha = 0.0125;
beta = 0.45;
T = 1000;

S0 = 20;
I0 = N - S0;

numberofrealisations=1;

for i=1:numberofrealisations
    S=S0;
    time=0;
    k=1;
    Aplot(1,i)=S;
    timeplot(1,i)=0;
    while (time < T)
        rr=rand(2,1);
        a0=alpha*S*(N-S)+beta*(N-S);
        time=time+(1/a0)*log(1/rr(1));
        if ( rr(2)*a0 < alpha*S*(N-S) )
            S=S-1;
        else
            S=S+1;
        end
        k=k+1;
        Aplot(k,i) = N-S; % calculate number infected
        timeplot(k,i)=time;
    end
end

% solve the deterministic (law of mass action) model
[tdet,Adet] = ode45(@(t,z) z*(alpha*N-beta-alpha*z),[0 T], I0);

figure(5);
for i = 1:numberofrealisations
    h=stairs(timeplot(:,i),Aplot(:,i)/N);
    set(h,'Linewidth',2);
    hold on;
end
plot(tdet,Adet/N,'--k','Linewidth',4);
ylabel('proportion of population infected');
xlabel('time');
xlim([0 50]);
ylim([0 1]);
set(gca,'Fontsize',20);
grid on;

str = ['R_0 = ',num2str(N*alpha/beta)];
title(str);
