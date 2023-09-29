function Figure1_5_a_b

% Matlab code for Figures 1.5(a) and 1.5(b) in book:
% Radek Erban and Jon Chapman, Stochastic Modelling of Reaction-Diffusion
% Processes, Cambridge University Press, 2019
 
close all;
rand('state',8);

k1=0.001;
k2=0.01;
k3=1.2;
k4=1;

Ainitial=0;
Binitial=0;
numberofrealisations=5;

for i=1:numberofrealisations

    A=Ainitial;
    B=Binitial;
    time=0;
    k=1;
    Aplot(1,i)=A;
    Bplot(1,i)=B;
    timeplot(1,i)=0;
    
    while (k<10000)
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

[tdet,det] = ode45(@myode,[0 1000],[Ainitial;Binitial]);

figure(1);
set(gca,'Fontsize',18);
plot(tdet,det(:,1),'--k','Linewidth',4);
hold on;
h=stairs(timeplot(:,1),Aplot(:,1));
set(h,'Color','r','Linewidth',1);
h=stairs(timeplot(:,2),Aplot(:,2));
set(h,'Color','m','Linewidth',1);
h=stairs(timeplot(:,3),Aplot(:,3));
set(h,'Color','b','Linewidth',1);
h=stairs(timeplot(:,4),Aplot(:,4));
set(h,'Color','g','Linewidth',1);
h=stairs(timeplot(:,5),Aplot(:,5));
set(h,'Color','c','Linewidth',1);
plot(tdet,det(:,1),'--k','Linewidth',4);
xlabel('time [sec]','interpreter','latex');
ylabel('number of $A$ molecules','interpreter','latex');
hh=legend('solution of ODEs');
set(hh,'interpreter','latex','Fontsize',18);
axis([0 100 0 25]);
set(gca,'Fontsize',18);

figure(2);
set(gca,'Fontsize',18);
plot(tdet,det(:,2),'--k','Linewidth',4);
hold on;
h=stairs(timeplot(:,1),Bplot(:,1));
set(h,'Color','r','Linewidth',1);
h=stairs(timeplot(:,2),Bplot(:,2));
set(h,'Color','m','Linewidth',1);
h=stairs(timeplot(:,3),Bplot(:,3));
set(h,'Color','b','Linewidth',1);
h=stairs(timeplot(:,4),Bplot(:,4));
set(h,'Color','g','Linewidth',1);
h=stairs(timeplot(:,5),Bplot(:,5));
set(h,'Color','c','Linewidth',1);
plot(tdet,det(:,2),'--k','Linewidth',4);
xlabel('time [sec]','interpreter','latex');
ylabel('number of $B$ molecules','interpreter','latex');
hh=legend('solution of ODEs');
set(hh,'interpreter','latex','location','northwest','Fontsize',18);
axis([0 100 0 25]);
set(gca,'Fontsize',18);

function dydt = myode(t,z)
k1=0.001;
k2=0.01;
k3=1.2;
k4=1;
dydt = [k3-2*k1*z(1)*z(1)-k2*z(1)*z(2); k4-k2*z(1)*z(2)];

