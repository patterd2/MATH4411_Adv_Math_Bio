function CH4_FP_reflect

dt=0.1;
finaltime=120*60;  % in minutes
D=0.0001;         % 10^{-4} mm^2 sec^{-1
numberofmolecules = 2500;
position=0.4*ones(numberofmolecules,1);
time=0;

while (time<finaltime)
       time=time+dt; 
       position=position+sqrt(2*D*dt)*randn(numberofmolecules,1);
       for i=1:numberofmolecules
           if (position(i)<0)
               position(i)=-position(i);
           end
           if (position(i)>=1)
               position(i)=2-position(i);
           end
       end
end

hist=zeros(40,1);
x=0.0125:0.025:0.9875;

for i=1:numberofmolecules
    hist(round(40*position(i)+0.5))=hist(round(40*position(i)+0.5))+1;
end

boxes=101;
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
axis([0 1 0 15+height]);
str = ['time = ', num2str(finaltime/60),' min'];
text(0.028,10+height,str,'Fontsize',22,'interpreter','tex');
hold on;
bar(x,hist,'b');
line(0:0.01:1,histdet,'Color','r','Linewidth',3);
ylabel('molecules per compartment');
xlabel('\eta_x');
legend('Fokker-Planck','SSA approx.');
set(gca,'Fontsize',20);
grid on;
box on;



