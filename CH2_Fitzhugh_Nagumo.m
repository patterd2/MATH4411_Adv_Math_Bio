%% Parameters and solution (from previous code)
a = 2;
I = 0.0;
c = -0.5;
epsilon = 0.1;
b = 0.5;
params = [a, b, c, I, epsilon];

V0 = 0.0;
W0 = -0.25;
y0 = [V0, W0];
tspan = [0 200];

[t, y] = ode45(@(t,y) fhn_model_epsilon(t, y, params), tspan, y0);

figure;
hold on;
plot(t, y(:,1), 'b-', t, y(:,2), 'r-', 'LineWidth', 2);
xlabel('Time');
ylabel('Variables');
legend('V (membrane potential)', 'W (recovery)',...
    'Location', 'northeast');
title('FitzHugh-Nagumo Model','FontWeight','normal');
grid on;
set(gca, 'FontSize', 20);

%% Phase Plane: Nullclines & Direction Field
V_min = -5; V_max = 5;
W_min = -5; W_max = 5;
Vvec = linspace(V_min,V_max,80);
Wvec = linspace(W_min,W_max,80);
[V, W] = meshgrid(Vvec,Wvec);

dVdt = V.*(V-a).*(1-V) - W + I;
dWdt = epsilon .* (V + c - b*W);
L = sqrt(dVdt.^2 + dWdt.^2);
dVdt_unit = dVdt ./ L;
dWdt_unit = dWdt ./ L;

% Nullclines
V_nullcline = Vvec.*(Vvec-a).*(1-Vvec) + I; 
W_nullcline = (Vvec+ c)/b;         

figure;
hold on;
quiver(V, W, dVdt_unit, dWdt_unit, 0.45, 'k'); % Direction field
plot(Vvec, V_nullcline, 'b', 'LineWidth', 3); % V-nullcline (dV/dt=0)
plot(Vvec, W_nullcline, 'r', 'LineWidth', 3);     % W-nullcline (dW/dt=0)
plot(y(:,1), y(:,2), 'g', 'LineWidth', 3);      % Trajectory
xlabel('V (membrane potential)');
ylabel('W (recovery variable)');
title('Fitzhugh-Nagumo nullclines and direction field',...
    'FontWeight','normal');
scatter(V0,W0,100,'k','filled');
legend({'Direction Field','V-nullcline','W-nullcline','Trajectory','Initial Condition'},...
    'Location','best');
set(gca,'FontSize',20);
xlim([-1 2.5])
ylim([-0.75 2.5])
grid on;

%% ODE Definition
function dydt = fhn_model_epsilon(~, y, p)
    V = y(1);
    W = y(2);
    a = p(1);
    b = p(2);
    I = p(4);
    c = p(3);
    epsilon = p(5);

    dVdt = V*(V-a)*(1-V) - W + I;
    dWdt = epsilon * (V + c - b*W);
    dydt = [dVdt; dWdt];
end
