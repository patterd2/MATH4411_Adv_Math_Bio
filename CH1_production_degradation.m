function CH1_production_degradation

k1=0.1;
k2=1;

Ainitial=100;
numberofrealisations=5;

for i=1:numberofrealisations
    A=Ainitial;
    time=0;
    k=1;
    Aplot(1,i)=A;
    timeplot(1,i)=0;
    while (k<1000)
          rr=rand(2,1);
          a0=k1*A+k2;
          time=time+(1/a0)*log(1/rr(1));
          if (rr(2)*a0<k1*A)
              A=A-1;
          else
              A=A+1;
          end
          k=k+1;
          Aplot(k,i)=A;
          timeplot(k,i)=time;
    end
end

tdet=0:0.2:100;
Adet(:)=(Ainitial-k2/k1)*exp(-k1*tdet(:))+k2/k1;

%% Plotting
figure;
set(gca,'Fontsize',18);
plot(tdet,Adet,'--k','Linewidth',4);
hold on
for j = 1:numberofrealisations
stairs(timeplot(:,j),Aplot(:,j),'Linewidth',2);
end
plot(tdet,Adet,'--k','Linewidth',4);
xlabel('time [sec]');
ylabel('number of molecules');
hh=legend('(deterministic) mean');
axis tight;
xlim([0 100]);
set(gca,'linewidth',1.5);
set(gca,'FontSize',20);
grid on;


