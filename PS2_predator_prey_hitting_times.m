function PS2_predator_prey
tic;
% set system parameters
b = 1;
p = 0.015;
d = 1;
T = 200; % length of time interval for each simulation

X0_list = 50; % initial conditions for system
Y0_list = 50;

%% Plot the level sets of the deterministic model for these parameters
% (thereby telling us what the solutions looks like in phase space)
V=@(x,y) p*x + p*y - d*log(x) - b*log(y);
[X,Y] = meshgrid(0:1:300);
z = V(X,Y);
figure(1);
scatter(X0_list,Y0_list,'r','LineWidth',5); % plot the initial condition
hold on;
contour(X,Y,z,-10:0.25:-4.5,'LineWidth',2);
set(gca,'Fontsize',20);
grid on;
title('Level sets of V(x,y)');
xlabel('x (prey)');
ylabel('y (prey)');
colormap("abyss");
legend('initial condition for sim');

%% Simulate the stochastic reaction directly using the Gillespie algorithm

numberofrealisations = 100; % increase for more accurate estimate of phi(n)

extinction_prey = 0;
extinction_preds = 0;

for i=1:numberofrealisations
    X = X0_list;
    Y = Y0_list;
    time=0;
    k=1;
    Xplot(1,i) = X;
    Yplot(1,i) = Y;
    timeplot(1,i)=0;
    while (time < T && X>0 && Y > 0 && max(X,Y) < 10^4)
        % simulations stops if both species go extinct (or become too
        % large)
        rr=rand(2,1);
        a0 = b*X;
        a1 = d*Y;
        a2 = p*X*Y;
        atot = a0+a1+a2; % propensity function
        time=time+(1/atot)*log(1/rr(1)); % next reaction time
        if ( rr(2)*atot < a0 )
            X = X + 1;
        elseif ( rr(2)*atot < a0 + a1 )
            Y = Y-1;
        else
            X = X - 1;
            Y = Y + 1;
        end
        k=k+1;
        Xplot(k,i) = X;
        Yplot(k,i) = Y;
        timeplot(k,i)=time;
        % disp(time);
        if X == 0
            extinction_prey(end+1) = time;
        end
        if Y == 0
            extinction_preds(end+1) = time;
        end
    end
end
extinction_preds = extinction_preds(2:end);
extinction_prey = extinction_prey(2:end);
%% Plot extinction time distributions
figure(2);
histogram(extinction_preds,floor(sqrt(length(extinction_preds))),'Normalization','pdf');
title('Predator Extinction Time Distribution');
set(gca,'Fontsize',20);
grid on;

figure(3);
histogram(extinction_prey(2:end),floor(sqrt(length(extinction_prey))),'Normalization','pdf');
title('Prey Extinction Time Distribution');
set(gca,'Fontsize',20);
grid on;

%% Plotting
for i = 1:min(3,numberofrealisations) % plot only the first 3 paths
    figure(4);
    [~, t2] = max(timeplot(:,i));
    stairs(timeplot(1:t2,i),Xplot(1:t2,i),'Linewidth',2);
    hold on;
    set(gca,'Fontsize',20);
    grid on;
    title('Prey Population');
    figure(5);
    [~, t2] = max(timeplot(:,i));
    stairs(timeplot(1:t2,i),Yplot(1:t2,i),'Linewidth',2);
    hold on;
    set(gca,'Fontsize',20);
    grid on;
    title('Predator Population');
end
toc;


