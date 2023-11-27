function CH3_CFP_HT1

k1=0.1;
k2=1;

A=k2/k1;
numberofrealisations=100000;
time=0;

p=zeros(22,1);

for i=1:numberofrealisations
    
    while (time<i)
          rr=rand(2,1);
          a0=k1*A+k2;
          time=time+(1/a0)*log(1/rr(1));
          if (rr(2)*a0<k1*A)
              A=A-1;
          else
              A=A+1;
          end
    end
    if ((A>0)&(A<23))
       p(A)=p(A)+1;
    end
end

p=p/numberofrealisations;

av=0:0.25:25;
pFP=exp(-2*av+((4*k2/k1)-1)*log(k1*av+k2));
pFP=pFP/(0.25*sum(pFP));

figure;
bar(p);
hold on;
plot(av,pFP,'r','Linewidth',3);
xlabel('number of molecules');
ylabel('stationary distribution');
legend('Gillespie SSA','Fokker-Planck');
set(gca,'XTick',[0 2 4 6 8 10 12 14 16 18 20 22]);
axis tight;
set(gca,'Fontsize',20);
grid on;

