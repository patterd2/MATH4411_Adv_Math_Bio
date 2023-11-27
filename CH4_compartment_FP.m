function CH4_compartment_FP

D=0.0001;         % 10^{-4} mm^2 sec^{-1}
finaltime = 1*60;   % in minutes
numberofmolecules=1500;

n=40; % number of compartments
h=1/n;
mesh=h/2:h:1-h/2;

A=zeros(n,1);

A(17)=numberofmolecules;

time=0;

while (time<finaltime)
       rr=rand(2,1);
       a0=2*D/(h*h)*(numberofmolecules-A(1)/2-A(n)/2);          
       time=time+(1/a0)*log(1/rr(1));
       ss=0;
       k=0;
       while ((ss<=rr(2)*a0)&(k<n-1))
              k=k+1;
              ss=ss+D/(h*h)*A(k);
       end
       if (ss>rr(2)*a0) 
           A(k)=A(k)-1;
           A(k+1)=A(k+1)+1;
       else
           k=1;
           while ((ss<=rr(2)*a0)&(k<n))
                  k=k+1;
                  ss=ss+D/(h*h)*A(k);
           end
           A(k)=A(k)-1;
           A(k-1)=A(k-1)+1;
       end
end

boxes=101;
dt=0.1;
histdet=zeros(boxes,1);
histdet(41)=numberofmolecules*100/40;
hh=dt*D*(boxes-1)*(boxes-1)/1/1;
for t=0:dt:finaltime
    histdetold=histdet;
    ii=1;
       histdet(ii)=histdetold(ii)+hh*(histdetold(ii+1)-histdetold(ii));
    for ii=2:boxes-1
       histdet(ii)=histdetold(ii)+hh*(histdetold(ii+1)+histdetold(ii-1)-2*histdetold(ii));
    end
    ii=boxes;
       histdet(ii)=histdetold(ii)+hh*(histdetold(ii-1)-histdetold(ii));
end

figure;
line(0:0.01:1,histdet,'Color','r','Linewidth',3);
height = ceil(max(histdet));
axis([0 1 0 25+height]);
str = ['time = ', num2str(finaltime/60),' min'];
text(0.028,20+height,str,'Fontsize',22,'interpreter','tex');
hold on;
bar(mesh,A,'b');
line(0:0.01:1,histdet,'Color','r','Linewidth',3);
xlabel('\eta_x','interpreter','tex');
ylabel('molecules per compartment','interpreter','tex');
legend('Fokker-Planck','Compartmental SSA');
set(gca,'Fontsize',20);
grid on;
box on;



