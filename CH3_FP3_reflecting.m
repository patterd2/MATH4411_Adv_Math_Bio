%function CH3_FP3_reflecting
tic;
%% parameters
dt=0.001;
finaltime=200;

n=finaltime/dt + 1;
numberofrealizations = 500;
L = 1;

%% time stepping to simualate the processes
xi = randn(n,numberofrealizations);
X = zeros(n,numberofrealizations);
Y = zeros(n,numberofrealizations);
tau = finaltime*ones(1,numberofrealizations);
f = @(x) x.^3;%-x.*(x-0.5).*(x-1);
g = @(x) 1/(2*sqrt(2));%-x.*(x-0.5).*(x-1);

for j = 1:numberofrealizations
    for i = 2:n
        %X(i,j) = X(i-1,j) + f(X(i-1,j))*dt + g(X(i-1,j))*sqrt(dt)*xi(i,j);
        Y(i,j) = Y(i-1,j) + f(Y(i-1,j))*dt + g(Y(i-1,j))*sqrt(dt)*xi(i,j);
        % if Y(i,j) < 0
        %    Y(i,j) = - Y(i-1,j) + Y(i-1,j)*dt - sqrt(2)*sqrt(dt)*xi(i,j);
        % end
        if Y(i,j) > L
            Y(i,j) = -2 + Y(i,j);
           %Y(i,j) = Y(i-1,j) + 2*(L-Y(i-1,j)) - f(Y(i-1,j))*dt - g(Y(i-1,j))*sqrt(dt)*xi(i,j);
        end
        if Y(i,j) < -L
           Y(i,j) = 2 + Y(i,j);
            %Y(i,j) = - 2*L - Y(i-1,j) - f(Y(i-1,j))*dt - g(Y(i-1,j))*sqrt(dt)*xi(i,j);
        end
    end
end
%% Plotting
bins = -2:0.05:2;
xs = -2:0.01:2;
xs2 = 0:0.01:L;
figure(1);
histogram(Y(end,:),bins,'Normalization','pdf');
hold on;
histogram(X(end,:),bins,'Normalization','pdf');
%plot(xs,exp(-0.5*xs.^2)./sqrt(2*pi),'-.g','LineWidth',4);
%plot(xs2,exp(-0.5*xs2.^2)/(sqrt(2*pi)*(normcdf(L)-0.5)),'-.m','LineWidth',4);
xlabel('x');
ylabel('p(x,100) [estimated]');
legend('reflecting BC','no BC');
set(gca,'Fontsize',20);
grid on;

figure(2);
%plot(0:dt:finaltime,Y(:,1),'-b','LineWidth',2);
scatter(0:dt:finaltime,Y(:,1));
hold on;
%plot(0:dt:finaltime,X(:,1),'-r','LineWidth',2);
scatter(0:dt:finaltime,X(:,1));
yline(0,'-.k','LineWidth',2);
xlabel('time');
legend('X_R(t) (reflected)','X(t)');
set(gca,'Fontsize',20);
grid on;

%%
toc;
