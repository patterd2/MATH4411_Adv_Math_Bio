%function CH2_no_stochastic_focusing

k1=100;
k2=1000;    
k3=0.01;  
k4=k2*99/10;   

for i=1:5

  A=10;
  B=100;
  C=0;
  time=0;
    
  for kk=1:2100
      while (time<kk*1)
           if (time>600) 
               A=5;
           end   
           rr=rand(2,1);
           a0=k1+k2*C+k3*B+k4*A*C;
           tau=(1/a0)*log(1/rr(1));
           ss=k1;     
           if (ss>rr(2)*a0)
              C=C+1;
           else
              ss=ss+k2*C;
              if (ss>rr(2)*a0)
                 C=C-1;
                 B=B+1;
              else
                 ss=ss+k3*B;
                 if (ss>rr(2)*a0)
                    B=B-1;
                 else
                    ss=ss+k4*A*C;
                    if (ss>rr(2)*a0)
                       C=C-1;
                    end
                 end
              end
           end
           time=time+tau;
      end
      AA(kk,i)=A;
      BB(kk,i)=B;
      CC(kk,i)=C;
      tt(kk,i)=time/60;
  end
end

[t,z] = ode15s(@myode,[0 2100],[100; 0]);

t=t/60;

figure(5);
set(gca,'Fontsize',18);
line([0 10 10 35],[10 10 5 5],'Color','b','Linewidth',4);
box on;
xlabel('time [min]');
ylabel('number of A molecules');
axis([0 35 0 22]);
set(gca,'Fontsize',20);
grid on;

figure(6);
plot(t,z(:,1),'k','Linewidth',4);
hold on;
for j = 1:5
    [~, t2] = max(tt(:,j));
    stairs(tt(1:t2,j),BB(1:t2,j),'Linewidth',2);
end
plot(t,z(:,1),'k','Linewidth',4);
xlabel('time [min]');
ylabel('number of B molecules');
hh=legend('solution of ODEs');
axis([0 35 50 400]);
set(gca,'Fontsize',20);
grid on;

function dydt = myode(t,z)
k1=100;
k2=1000;
k3=0.01;
k4=k2*99/10;
A=10;
if (t>600)
   A=5;
end
dydt =[k2*z(2)-k3*z(1); k1-k2*z(2)-k4*A*z(2)];
end
