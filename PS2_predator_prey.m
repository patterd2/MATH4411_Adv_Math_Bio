function PS2_predator_prey

% set system parameters
b = 1;
p = 0.015;
d = 1;
T = 100; % length of time interval for each simulation

X0_list = 50; % initial conditions for system
Y0_list = 50;

%% Simulate the stochastic reaction directly using the Gillespie algorithm

numberofrealisations = 1;

for i=1:numberofrealisations
    X = X0_list;
    Y = Y0_list;
    time=0;
    k=1;
    Xplot(1,i) = X;
    Yplot(1,i) = Y;
    timeplot(1,i)=0;
    while (time < T && (X>0 || Y >0) && max(X,Y) < 10^5) 
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
    end
end

%% Plotting
for i = 1:numberofrealisations
    figure(1);
    [~, t2] = max(timeplot(:,i));
    stairs(timeplot(1:t2,i),Xplot(1:t2,i),'Linewidth',2);
    hold on;
    set(gca,'Fontsize',20);
    grid on;
    title('Prey Population');
    figure(2);
    [~, t2] = max(timeplot(:,i));
    stairs(timeplot(1:t2,i),Yplot(1:t2,i),'Linewidth',2);
    hold on;
    set(gca,'Fontsize',20);
    grid on;
    title('Predator Population');
end



