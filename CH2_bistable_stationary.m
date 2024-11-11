function CH2_bistable_stationary

k1=0.00025;
k2=0.18;
k3=36;
k4=2200;

X=0;
time=0;
phist=zeros(601,1);

for i=1:10000   
    while (time<i*0.1)
          rr=rand(2,1);
          a0=k1*X*(X-1)*(X-2)+k2*X*(X-1)+k3*X+k4;
          tau=(1/a0)*log(1/rr(1));
          if (rr(2)*a0<(k2*X*(X-1)+k4))
             X=X+1;
          else
             X=X-1;
          end  
          time=time+tau;
    end
    if (X<601.5)
       phist(X+1)=phist(X+1)+1;
    end
end

phist=phist/sum(phist);

pCME(1)=k4/k3;
f(1)=k4;
g(1)=k3;
sumpCME=1+pCME(1);

for ii=2:600
    f(ii)=k2*ii*(ii-1)+k4;
    g(ii)=k1*ii*(ii-1)*(ii-2)+k3*ii;
    pCME(ii)=pCME(ii-1)*f(ii-1)/g(ii);
    sumpCME=sumpCME+pCME(ii);
end

pCME=pCME/sumpCME;

figure;
bar([700],[0],'b');
hold on;
plot(pCME,'r','Linewidth',3);
stem([0:1:600],phist,'b','Markersize',0.02);
plot(pCME,'r','Linewidth',3);
set(gca,'XTick',[0 100 200 300 400 500 600]);
box on;
xlabel('number of molecules');
legend('Gillespie SSA','Stationary Distribution');
set(gca,'Fontsize',20);
grid on;
axis tight;
xlim([0 600]);
