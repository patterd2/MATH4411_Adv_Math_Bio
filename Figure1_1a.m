function Figure1_1a

% Matlab code for Figure 1.1(a) in book:
% Radek Erban and Jon Chapman, Stochastic Modelling of Reaction-Diffusion
% Processes, Cambridge University Press, 2019

close all;
rand('state',22);

k1=0.1;
Ainitial=20;
dt=0.005;
numberofrealisations=10;

for i=1:numberofrealisations

    A=Ainitial;
    time=0;
    k=1;
    Aplot(1,i)=A;
    timeplot(1,i)=0;

    while (k<Ainitial+1)
         rr=rand(1,1);
         time=time+dt;
         if (rr<k1*A*dt)
             A=A-1;
             k=k+1;
             Aplot(k,i)=A;
             timeplot(k,i)=time;
         end
    end
    Aplot(k+1,i)=A;
    if (time<30)
        timeplot(k+1,i)=30;
    else
        timeplot(k+1,i)=time;
    end
end

tdet=zeros(301,1);
Adet=zeros(301,1);
tdet=[0:0.1:30];
Adet(:)=Ainitial*exp(-k1*tdet(:));

figure(1);
set(gca,'Fontsize',18);
h=stairs(timeplot(:,1),Aplot(:,1));
set(h,'Color','b','Linewidth',2);
hold on;
h=stairs(timeplot(:,2),Aplot(:,2));
set(h,'Color','r','Linewidth',2);
xlabel('time [sec]','interpreter','latex');
ylabel('number of molecules','interpreter','latex');
hh=legend('first realization','second realization');
set(hh,'interpreter','latex','Fontsize',18);
axis([0 30 0 (Ainitial+1)]);
set(gca,'Fontsize',18);


