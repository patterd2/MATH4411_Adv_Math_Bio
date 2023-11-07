function CH3_example_3

dt=0.001;
T=50; 
numberofrealisations=5;


n=T/dt;
X=zeros(numberofrealisations,n);
time=0:dt:T;

for i=1:numberofrealisations
    dX = dt*ones(1,n)+sqrt(dt)*randn(1,n);   
    X(i,:)=cumsum(dX(1,:));             
end

%% Plotting
figure(1);
for j=1:numberofrealisations
plot(time,[0,X(j,:)]);   
hold on;
end
line([0 T],[0 T],'Linewidth',2,'Linestyle','--','Color','k');
xlabel('t','interpreter','latex');
ylabel('X(t)');
set(get(gca,'ylabel'),'rotation',0);
axis tight;
set(gca,'Fontsize',20);
grid on;

