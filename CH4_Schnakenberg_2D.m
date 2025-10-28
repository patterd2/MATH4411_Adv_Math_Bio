% CH4_Schnakenberg_2D.m
% Schnakenberg stochastic reaction-diffusion model on a 2D spatial domain

%% Parameters
DA = 1.0e-5;
DB = 1.0e-3;
Kx = 40;              % grid size in x-direction
Ky = 40;              % grid size in y-direction
Lx = 1;               % domain size x
Ly = 1;               % domain size y
T = 2000;             % simulation end time
M = 2e7;              % max number of events

h = Lx / Kx;          % assuming square grid
k1 = 1.0e-6;
k2 = 1;
k3 = 0.02;
k4 = 3;

%% Initialization
A = round(200 * ones(Kx, Ky));
B = round(70 * ones(Kx, Ky));
time = 0;
kk = 1;

%% Diffusion rates per molecule
dA = DA / h^2;
dB = DB / h^2;

%% SSA simulation
while time < T && kk < M
    Ak = A;
    Bk = B;

    % Reaction propensities
    alpha1 = k1 * Ak .* max(Ak - 1, 0) .* Bk;
    alpha2 = k2 * ones(Kx, Ky);
    alpha3 = k3 * Ak;
    alpha4 = k4 * ones(Kx, Ky);

    % Diffusion propensities per site
    diffA = dA * Ak;
    diffB = dB * Bk;

    % Total propensities
    alphas = [sum(alpha1(:)), sum(alpha2(:)), sum(alpha3(:)), sum(alpha4(:)), sum(diffA(:)), sum(diffB(:))];
    total_alpha = sum(alphas);
    tau = (1 / total_alpha) * log(1 / rand);

    % Choose event
    r = rand * total_alpha;
    thresholds = cumsum(alphas);
    event = find(r <= thresholds, 1);

    % Update system
    switch event
        case 1 % A + A + B -> 3A
            idx = find(alpha1 > 0);
            chosen = idx(randi(length(idx)));
            [i,j] = ind2sub([Kx, Ky], chosen);
            if Ak(i,j) >= 2 && Bk(i,j) >= 1
                A(i,j) = Ak(i,j) + 1;
                B(i,j) = Bk(i,j) - 1;
            end
        case 2 % feed A
            [i,j] = ind2sub([Kx, Ky], randi(Kx * Ky));
            A(i,j) = Ak(i,j) + 1;
        case 3 % A -> nothing
            idx = find(alpha3 > 0);
            chosen = idx(randi(length(idx)));
            [i,j] = ind2sub([Kx, Ky], chosen);
            if Ak(i,j) >= 1
                A(i,j) = Ak(i,j) - 1;
            end
        case 4 % feed B
            [i,j] = ind2sub([Kx, Ky], randi(Kx * Ky));
            B(i,j) = Bk(i,j) + 1;
        case 5 % diffuse A
            idx = find(diffA > 0);
            chosen = idx(randi(length(idx)));
            [i,j] = ind2sub([Kx, Ky], chosen);
            [ni,nj] = random_neighbor(i, j, Kx, Ky);
            A(i,j) = Ak(i,j) - 1;
            A(ni,nj) = Ak(ni,nj) + 1;
        case 6 % diffuse B
            idx = find(diffB > 0);
            chosen = idx(randi(length(idx)));
            [i,j] = ind2sub([Kx, Ky], chosen);
            [ni,nj] = random_neighbor(i, j, Kx, Ky);
            B(i,j) = Bk(i,j) - 1;
            B(ni,nj) = Bk(ni,nj) + 1;
    end

    time = time + tau;
    kk = kk + 1;

    if mod(kk,50000) == 0
        % Plot solution at each event
        figure(1); clf;
        subplot(1,2,1);
        imagesc(A, [100 400]);
        title(['A concentration, time = ', num2str(time, '%.2f')]);
        axis square; colorbar;

        subplot(1,2,2);
        imagesc(B, [0 120]);
        title('B concentration');
        axis square; colorbar;
        drawnow;
    end
end

%% Final snapshot
figure;
subplot(1,2,1);
imagesc(A, [0 500]);
title('A concentration');
axis square; colorbar;

subplot(1,2,2);
imagesc(B, [0 120]);
title('B concentration');
axis square; colorbar;

%% Helper function
function [nx, ny] = random_neighbor(x, y, Kx, Ky)
dirs = [0 1; 1 0; 0 -1; -1 0]; % right, down, left, up
dir = dirs(randi(4), :);
nx = x + dir(1);
ny = y + dir(2);
if nx < 1, nx = 1; elseif nx > Kx, nx = Kx; end
if ny < 1, ny = 1; elseif ny > Ky, ny = Ky; end
end
