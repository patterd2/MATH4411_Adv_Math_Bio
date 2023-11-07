function CH3_example_4

k1=0.001;
k2=0.75;
k3=165;
k4=10000;
k5=200;

dt=0.01;
finaltime=800; 
n=finaltime/dt;

X=zeros(n+1,1);
time=0:dt:finaltime;

for i=1:n
    X(i+1) = X(i)+dt*(-k1*X(i)*X(i)*X(i)+k2*X(i)*X(i)-k3*X(i)+k4)+k5*sqrt(dt)*randn(1,1);   
end

[t,z] = ode15s(@(t,z) -k1*z(1)*z(1)*z(1)+k2*z(1)*z(1)-k3*z(1)+k4,[0 finaltime],0);

figure(1);
plot(time,X,'b');   
hold on;
plot(t,z,'r','Linewidth',3);
xlabel('t');
ylabel('X(t)');
set(get(gca,'ylabel'),'rotation',0);
legend('SDE','ODE');
axis tight
box on;
set(gca,'Fontsize',20);
grid on;
