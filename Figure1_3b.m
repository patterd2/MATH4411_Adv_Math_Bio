function Figure1_3b

% Matlab code for Figure 1.3(b) in book:
% Radek Erban and Jon Chapman, Stochastic Modelling of Reaction-Diffusion
% Processes, Cambridge University Press, 2019

close all;
rand('state',10);

k1=0.005;
k2=1;

A=sqrt(k2/(2*k1));
numberofrealisations=100000;
time=0;

p=zeros(24,1);

for i=1:numberofrealisations
    
    while (time<i)
          rr=rand(2,1);
          a0=k1*A*(A-1)+k2;
          time=time+(1/a0)*log(1/rr(1));
          if (rr(2)*a0<k1*A*(A-1))
              A=A-2;
          else
              A=A+1;
          end
    end
    if (A<23.5)
       p(A+1)=p(A+1)+1;
    end
end

p=p/numberofrealisations;

% p(n+1) is the stationary probability that A=n for n=0,1,2,3,.....

for n=0:25
    phi(n+1)=((sqrt(k2/k1)^n)/factorial(n))*besseli(n-1,2*(sqrt(k2/k1)));
end
phi=phi/sum(phi);

figure(1);
set(gca,'Fontsize',18);
h1=bar([0:23],p);
set(h1,'FaceColor',[0.88 0.88 0.88]);
hold on;
h2=plot([0:25],phi,'r','Linewidth',3);
xlabel('number of molecules','interpreter','latex');
ylabel('stationary distribution','interpreter','latex');
hh=legend([h1 h2],'Gillespie SSA','analysis');
set(hh,'interpreter','latex','Fontsize',18);
set(gca,'XTick',[0 5 10 15 20]);
axis([0 22.5 0 1.05*max(p)]);
set(gca,'Fontsize',18);
