function CH3_FP1

dt=0.001;
finaltime=1; 

n=finaltime/dt;
numberofrealizations=1000;
xg=-5:0.05:5;

x=-4.9:0.2:4.9;
hist=zeros(length(x),1);

for j=1:10
    dX = sqrt(dt)*randn(n,numberofrealizations);   
    position=sum(dX);             
    for i=1:numberofrealizations
        hist(round(5*position(i)+25.5))=hist(round(5*position(i)+25.5))+1;
    end
end

hist=hist/numberofrealizations/2;

gaussian=exp(-xg.*xg/2)/(sqrt(2*pi));

%% Plotting
figure(1);
bar(x,hist);
hold on;
plot(xg,gaussian,'r','Linewidth',3);
xlabel('x');
ylabel('p(x,1)');
legend('SSA','Fokker-Planck');
set(gca,'XTick',[-4 -3 -2 -1 0 1 2 3 4]);
axis([-4 4 0 0.45]);
set(gca,'Fontsize',20);
grid on;
