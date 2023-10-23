function CH1_more_efficient_SSA_degradation

k=0.1; % reaction rate
Ainitial=20; % initial number of molecules
numberofrealisations=200; % number of sample paths
T = 30; % end time for simulation 

%% More efficient SSA
for i=1:numberofrealisations
    A=Ainitial;
    time=0;
    j=1;
    Aplot(1,i)=A;
    timeplot(1,i)=0;

    while (j<Ainitial+1)
         rr=rand(2,1);
         a0=k*A;
         time=time+(1/a0)*log(1/rr(1));
         A=A-1;
         j=j+1;
         Aplot(j,i)=A;
         timeplot(j,i)=time;
    end
    Aplot(j+1,i)=A;
    timeplot(j+1,i)=30;
end
%% Plotting of sample paths vs deterministic mean
tdet=0:0.1:T;
Adet(:)=Ainitial*exp(-k*tdet(:)); % compute deterministic mean for comparison

figure;
plot(tdet,Adet,'k--','Linewidth',4);
hold on
for i=1:numberofrealisations
    h=stairs(timeplot(:,i),Aplot(:,i));
    set(h,'Linewidth',2);
end
hold on;
plot(tdet,Adet,'k--','Linewidth',4);
xlabel('time [sec]');
ylabel('number of molecules');
hh=legend('deterministic mean');
axis([0 30 0 (Ainitial+1)]);
set(gca,'linewidth',1.5);
set(gca,'FontSize',20);
grid on;
%% Finite-time distributions
% Goal: Estimate the distribution of the process at time t = 10
% This part of the code needs a larger number of realisations to work well
FTvalues = zeros(1,numberofrealisations);
for i = 1:numberofrealisations
    [temp1, temp2] = max(timeplot(timeplot(:,i)<10,i)); % find the last jump before t = 10
    FTvalues(i) = Aplot(temp2,i); % extract the value at t = 10
end
figure;
edges = 0.5:1:20;
histogram(FTvalues,edges,'Normalization','pdf');
hold on;
xline(Ainitial*exp(-k*10),'k--','Linewidth',4)
xlim([0 Ainitial]);
legend('Estimated PMF (t=10)','det. mean (t = 10)');
set(gca,'linewidth',1.5);
set(gca,'FontSize',20);
grid on;