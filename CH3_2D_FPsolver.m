% Parameters
Lx = 1;                % Domain length in x
Ly = 1;                % Domain length in y
Nx = 25;               % Number of spatial points in x
Ny = 25;               % Number of spatial points in y
dx = Lx / (Nx - 1);    % Spatial step size in x
dy = Ly / (Ny - 1);    % Spatial step size in y
dt = 0.001;             % Time step size
Nt = 1000;             % Number of time steps

x = linspace(0, Lx, Nx)';
y = linspace(0, Ly, Ny)';
[X, Y] = meshgrid(x, y);

% Drift and diffusion terms (example: constant)
Ax = @(X,Y) 0;          % Drift in x
Ay = @(X,Y) 0;          % Drift in y
Dx = @(X,Y) 1;          % Diffusion in x
Dy = @(X,Y) 1;          % Diffusion in y

% Initial condition (Gaussian centered at midpoints)
P = exp(-50 * ((X - Lx/2).^2 + (Y - Ly/2).^2));
P = P / sum(P(:)) / (dx * dy); % Normalize initial condition

% Precompute coefficients
rx = dt / (2 * dx^2);
ry = dt / (2 * dy^2);

% Create sparse matrices for x and y directions
Ix = speye(Nx); Iy = speye(Ny);
DxMat = spdiags([-rx * ones(Nx,1), (1 + 2 * rx) * ones(Nx,1), -rx * ones(Nx,1)], [-1 0 1], Nx, Nx);
DyMat = spdiags([-ry * ones(Ny,1), (1 + 2 * ry) * ones(Ny,1), -ry * ones(Ny,1)], [-1 0 1], Ny, Ny);

% No-flux boundary conditions for DxMat and DyMat
DxMat(1,1) = 1; DxMat(1,2) = -1;
DxMat(end,end) = 1; DxMat(end,end-1) = -1;

DyMat(1,1) = 1; DyMat(1,2) = -1;
DyMat(end,end) = 1; DyMat(end,end-1) = -1;

% Time-stepping loop
for n = 1:Nt
    % Implicit solution in x, then in y
    for j = 1:Ny  % Solve for each row (x-direction)
        P(j,:) = (DxMat \ (DxMat * P(j,:)'))';
    end
    for i = 1:Nx  % Solve for each column (y-direction)
        P(:,i) = (DyMat \ (DyMat * P(:,i)));
    end
    % Plot the solution
    figure(1);
    hold on;
    surf(X, Y, P, 'EdgeColor', 'none');
    xlabel('x');
    ylabel('y');
    zlabel('P(x,y,t)');

    % Normalize
    P = P / sum(P(:)) / (dx * dy);

end



