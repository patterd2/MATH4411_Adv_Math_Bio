function CH3_FP3_reflecting
%% parameters
dt=0.001;
finaltime=10;

n=finaltime/dt + 1;
numberofrealizations=10000;

%% time stepping to simualate the processes
xi = randn(n,numberofrealizations);
X = zeros(n,numberofrealizations);
Y = zeros(n,numberofrealizations);

for j = 1:numberofrealizations
    for i = 2:n
        X(i,j) = X(i-1,j) - X(i-1,j)*dt + sqrt(dt)*xi(i,j);
        Y(i,j) = Y(i-1,j) - Y(i-1,j)*dt + sqrt(dt)*xi(i,j);
        if Y(i,j) < 0
            Y(i,j) = -Y(i-1,j) + Y(i-1,j)*dt - sqrt(dt)*xi(i,j);
        end
    end
end

%% Plotting
bins = -2.2:0.1:2.2;
figure(1);
histogram(Y(end,:),bins,'Normalization','probability');
hold on;
histogram(X(end,:),bins,'Normalization','probability');
xlabel('x');
ylabel('p(x,10) [estimated]');
legend('reflecting BC','no BC');
set(gca,'Fontsize',20);
grid on;

figure(2);
plot(0:dt:finaltime,Y(:,1),'-b','LineWidth',2);
hold on;
plot(0:dt:finaltime,X(:,1),'-r','LineWidth',2);
yline(0,'-.k','LineWidth',2);
xlabel('time');
legend('X_R(t) (reflected)','X(t)');
set(gca,'Fontsize',20);
grid on;
