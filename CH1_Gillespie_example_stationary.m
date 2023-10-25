function CH1_Gillespie_example_stationary

k1=0.001;
k2=0.01;
k3=1.2;
k4=1;

Ainitial=sqrt((k3-k4)/(2*k1));
Binitial=k4/(k2*Ainitial);
numberofrealisations=100000;
pa=zeros(101,1);
p2D=zeros(31,31);

A=round(Ainitial);
B=round(Binitial);
time=0;

for i=1:numberofrealisations
    while (time<i)
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
    end
    if (A<100.5)
       pa(A+1)=pa(A+1)+1;
    end
    if (A<30.5)
       if (B<30.5)
           p2D(A+1,B+1)=p2D(A+1,B+1)+1;
       end
    end
end

pa=pa/numberofrealisations;
p2D=p2D/numberofrealisations;
 
figure(1);
imagesc(0:30,0:30,p2D');
set(gca,'YDir','normal');
xlabel('number of A molecules');
ylabel('number of B molecules');
colorbar;
title('Estimated joint PMF');
axis([0 30 0 30]);
view(0,90);
set(gca,'Fontsize',20);

figure(2);
bar(0:100,pa);
xlabel('number of A molecules');
title('Stationary Distribution (estimated)');
set(gca,'XTick',[0 5 10 15 20]);
axis([-0.5 23.5 0 1.05*max(pa)]);
set(gca,'Fontsize',20);
grid on;

