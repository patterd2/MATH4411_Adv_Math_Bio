%% Parameters

alpha = 1; % Recruitment/invasion G to F
f_0 = 0.1;
f_1 = 0.9;
t_2 = 0.4;
s_2 = 0.05;

params = [alpha, f_0, f_1, t_2, s_2];

%% Compute the equilibrium solutions
f = @(x) (1-x).*(phi(x, f_0, f_1, t_2, s_2)-alpha*x);
x_grid = 0:0.001:1;
fx = f(x_grid);
figure;
hold on;
plot(x_grid, fx,'LineWidth',3);
str = sprintf('$\\alpha = %.1f$', alpha);
annotation('textbox', [0.02, 0.85, 0.15, 0.1], 'String', str, ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white', 'FontSize',...
    26, 'Interpreter', 'latex');
xlabel('G','Interpreter','latex');
yline(0,'--r','LineWidth',3)
ylabel('f(G)','Interpreter','latex','Rotation',0);
title('$\frac{d}{dt}G = f(G) := (1-G)(\phi(G)-\alpha G)$',...
    'Interpreter','latex')
set(gca,'FontSize',25);
grid on;
ylim([min(fx)-0.025 max(fx)+0.025]);

%% Solve ODE
% Initial conditions: G = 0.6, F = 1-G (cover fractions)
G0 = 0.2;
y0 = [G0, 1-G0];

% Time span
tspan = [0 100];

%%
% ODE solver
[t, y] = ode45(@(t, y) sl_model(t, y, params), tspan, y0);

%% Plotting
figure; 
plot(t, y, 'LineWidth', 4);
xlabel('Time'); ylabel('Cover Fraction');
legend('Grass (G)','Forest (F)','Location','best');
title('Savanna-Forest (ODE) Model');
grid on;
str = sprintf('$\\alpha = %.1f$', alpha);
annotation('textbox', [0.02, 0.85, 0.15, 0.1], 'String', str, ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white', 'FontSize',...
    26, 'Interpreter', 'latex');
set(gca,'FontSize',25);
ylim([-0.005 1.005]);

%% Phase space plot
% Phase space domain
% g_range = linspace(0, 1, 30);
% f_range = linspace(0, 1, 30);
% [G, F] = meshgrid(g_range, f_range);
% 
% dG = zeros(size(G));
% dF = zeros(size(F));
% 
% for i=1:numel(G)
%     if G(i) + F(i) <= 1
%         derivs = sl_model(0, [G(i); F(i)], params);
%         dG(i) = derivs(1);
%         dF(i) = derivs(2);
%     else
%         dG(i) = 0;
%         dF(i) = 0;
%     end
% end
% 
% L = sqrt(dG.^2 + dS.^2);               % Compute lengths
% dG_unit = dG ./ L;                     % Normalize
% dF_unit = dF ./ L;
% % To avoid division by zero, set NaN or zero where magnitude is zero
% dG_unit(L == 0) = 0;
% dF_unit(L == 0) = 0;
% 
% % Plot direction field (vector field) in phase space
% figure;
% quiver(G, F, dG_unit, dF_unit, 'k')
% hold on
% xlabel('Grass Fraction (G)')
% ylabel('Forest Fraction (F)')
% title('Phase Space Direction Field')
% xlim([0 1]); ylim([0 1]);
% str = sprintf('$\\alpha = %.1f$', alpha);
% annotation('textbox', [0.435, 0.75, 0.25, 0.15], 'String', str, ...
%     'FitBoxToText', 'on', 'BackgroundColor', 'white', 'FontSize',...
%     26, 'Interpreter', 'latex');
% set(gca,'FontSize',25);

%% function definitions
function dydt = sl_model(~, y, p)
G = y(1); F = y(2);
alpha = p(1);
f_0 = p(2);
f_1 = p(3);
t_2 = p(4);
s_2 = p(5);

dG = phi(G, f_0, f_1, t_2, s_2).*F - alpha*G.*F;
dF = alpha*G.*F - phi(G, f_0, f_1, t_2, s_2).*F;
dydt = [dG; dF];
end

function y = phi(g, f_0, f_1, t_2, s_2)
    y = f_0 + (f_1-f_0)./(1 + exp(-(g-t_2)/s_2));
end

figure(1);