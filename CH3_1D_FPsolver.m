% Parameters
L = 1;                 % Domain length
Nx = 100;              % Number of spatial points
dx = L / (Nx - 1);     % Spatial step size
x = linspace(0, L, Nx)'; % Spatial grid
dt = 0.0001;             % Time step size
Nt = 1000;             % Number of time steps

%%

% Drift and diffusion terms (user-defined)
A = @(x) 0;            % Drift term (example: zero drift)
D = @(x) 1;            % Diffusion term (example: constant diffusion)

% Initial condition
P = exp(-100 * (x - L/2).^2); % Gaussian centered at L/2
P = P / sum(P) / dx;          % Normalize initial condition

% Crank-Nicholson matrices
r = dt / (2 * dx^2); % Stability parameter, must be < 1/2
main_diag = (1 + 2 * r) * ones(Nx, 1);    % Main diagonal
upper_diag = -r * ones(Nx, 1);            % Upper diagonal (size Nx)
lower_diag = -r * ones(Nx, 1);            % Lower diagonal (size Nx)

% Set boundary conditions: Neumann (no-flux)
upper_diag(end) = 0;   % No connection beyond the last point
lower_diag(1) = 0;     % No connection before the first point

% Construct sparse matrices
A_matrix = spdiags([lower_diag main_diag upper_diag], [-1 0 1], Nx, Nx);
B_matrix = spdiags([-lower_diag (1 - 2 * r) * ones(Nx, 1) -upper_diag], [-1 0 1], Nx, Nx);

% No-flux boundary conditions (Neumann)
A_matrix(1,1) = 1; A_matrix(1,2) = -1;
A_matrix(end,end) = 1; A_matrix(end,end-1) = -1;

B_matrix(1,1) = 1; B_matrix(1,2) = -1;
B_matrix(end,end) = 1; B_matrix(end,end-1) = -1;

% Time-stepping loop
for n = 1:Nt
    % Solve the linear system for Crank-Nicholson step
    P = A_matrix \ (B_matrix * P);

    % Ensure probability conservation
    P = P / sum(P) / dx;
    figure(1);
    hold on;
    plot(x, P, 'LineWidth', 2);
    xlabel('x');
    ylabel('P(x,t)');
    title('Solution of the Fokker-Planck Equation');
end

% Plot the solution
figure;
plot(x, P, 'LineWidth', 2);
xlabel('x');
ylabel('P(x,t)');
title('Solution of the Fokker-Planck Equation');
