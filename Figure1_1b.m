function Figure1_1b

% Matlab code for Figure 1.1(b) in book:
% Radek Erban and Jon Chapman, Stochastic Modelling of Reaction-Diffusion
% Processes, Cambridge University Press, 2019

close all;
rand('state',6);

k1=0.1;
Ainitial=20;
numberofrealisations=10;

for i=1:numberofrealisations

    A=Ainitial;
    time=0;
    k=1;
    Aplot(1,i)=A;
    timeplot(1,i)=0;

    while (k<Ainitial+1)
         rr=rand(2,1);
         a0=k1*A;
         time=time+(1/a0)*log(1/rr(1));
         A=A-1;
         k=k+1;
         Aplot(k,i)=A;
         timeplot(k,i)=time;
    end
    
    Aplot(k+1,i)=A;
    timeplot(k+1,i)=30;
end

tdet=zeros(301,1);
Adet=zeros(301,1);
tdet=[0:0.1:30];
Adet(:)=Ainitial*exp(-k1*tdet(:));

figure(1);
set(gca,'Fontsize',18);
plot(tdet,Adet,'k--','Linewidth',4);
hold on
h=stairs(timeplot(:,6),Aplot(:,6));
set(h,'Color',[0.9 0.9 0],'Linewidth',1);
h=stairs(timeplot(:,1),Aplot(:,1));
set(h,'Color',[0.3 0.6 0.9],'Linewidth',1);
h=stairs(timeplot(:,2),Aplot(:,2));
set(h,'Color',[0.6 0.9 0.3],'Linewidth',1);
h=stairs(timeplot(:,3),Aplot(:,3));
set(h,'Color',[0.9 0.3 0.6],'Linewidth',1);
h=stairs(timeplot(:,4),Aplot(:,4));
set(h,'Color','m','Linewidth',1);
h=stairs(timeplot(:,5),Aplot(:,5));
set(h,'Color','c','Linewidth',1);
h=stairs(timeplot(:,7),Aplot(:,7));
set(h,'Color','k','Linewidth',1);
h=stairs(timeplot(:,8),Aplot(:,8));
set(h,'Color','g','Linewidth',1);
h=stairs(timeplot(:,9),Aplot(:,9));
set(h,'Color','r','Linewidth',1);
h=stairs(timeplot(:,10),Aplot(:,10));
set(h,'Color','b','Linewidth',1);
plot(tdet,Adet,'k--','Linewidth',4);
xlabel('time [sec]','interpreter','latex');
ylabel('number of molecules','interpreter','latex');
hh=legend('mean');
set(hh,'interpreter','latex','Fontsize',18);
axis([0 30 0 (Ainitial+1)]);
set(gca,'Fontsize',18);

