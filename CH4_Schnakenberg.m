function CH4_Schnakenberg

% Compartmental model of stochastic reaction-diffusion with patterns
% Schnakenberg reaction kinetics with discrete diffusive movement

%% Set parameter values
DA = 1.000e-05;
DB = 1.0000e-03;
K = 40; % number of compartments
L = 1; % spatial domain = [0,L]
T = 2000; % final time for the simulations
M = 20000000; % max number of events

h = L/K; % spatial discretisation parameter
dA = 2*(DA/h^2)*ones(1,K);
dA(1,1) = DA/h^2;
dA(1,K) = DA/h^2;
dB = 2*(DB/h^2)*ones(1,K);
dB(1,1) = DB/h^2;
dB(1,K) = DB/h^2;
k1 = 1.0e-06; % already scaled by volume to give effective rates
k2 = 1;
k3 = 0.02;
k4 = 3;

% uniform initial conditions
A0=round(200*ones(M,K));
B0=round(70*ones(M,K));

time=0;
kk=0; % event count
A = A0;
B = B0;

Atot = sum(A0(1,:));
Btot = sum(B0(1,:));

figure(1);
fig=gcf;
fig.Units='normalized';
fig.OuterPosition=[0 0 1 1];
%% SSA
while time < T && kk < M
    kk=kk+1;
    alpha1 = k1*A(kk,:).*max((A(kk,:)-1),0).*B(kk,:);
    alpha2 = k2*K;
    alpha3 = k3*A(kk,:);
    alpha4 = k4*K;
    alpha5 = sum(dA.*A(kk,:));
    alpha6 = sum(dB.*B(kk,:));
    alphas = [sum(alpha1) alpha2 sum(alpha3) alpha4 alpha5 alpha6];
    alpha0 = sum(alphas);
    tau=(1/alpha0)*log(1/rand);
    alpha_norm = cumsum(alphas)/alpha0;
    event = find(alpha_norm > rand,1,'first'); % select the event to occur
    %keyboard;
    A(kk+1,:) = A(kk,:);
    B(kk+1,:) = B(kk,:);
    if event == 1
        subevent = find(cumsum(alpha1)/sum(alpha1) > rand,1,'first');
        A(kk+1,subevent) = A(kk,subevent) + 1;
        B(kk+1,subevent) = B(kk,subevent) - 1;
    elseif event == 2
        subevent = randi([1 K]);
        A(kk+1,subevent) = A(kk,subevent) + 1;
    elseif event == 3
        subevent = find(cumsum(alpha3)/sum(alpha3) > rand,1,'first');
        A(kk+1,subevent) = A(kk,subevent) - 1;
    elseif event == 4
        subevent = randi([1 K]);
        B(kk+1,subevent) = B(kk,subevent) + 1;
    elseif event == 5
        subevent = find(cumsum(dA.*A(kk,:))/sum(dA.*A(kk,:)) > rand,1,'first');
        if subevent == 1
            A(kk+1,1) = A(kk,1) - 1;
            A(kk+1,2) = A(kk,2) + 1;
        elseif subevent == K
            A(kk+1,K) = A(kk,K) - 1;
            A(kk+1,K-1) = A(kk,K-1) + 1;
        else
            A(kk+1,subevent) = A(kk,subevent) - 1;
            left_right = sign(rand-0.5);
            A(kk+1,subevent+left_right) = A(kk,subevent+left_right) + 1;
        end
    elseif event == 6
        subevent = find(cumsum(dB.*B(kk,:))/sum(dB.*B(kk,:)) > rand,1,'first');
        if subevent == 1
            B(kk+1,1) = B(kk,1) - 1;
            B(kk+1,2) = B(kk,2) + 1;
        elseif subevent == K
            B(kk+1,K) = B(kk,K) - 1;
            B(kk+1,K-1) = B(kk,K-1) + 1;
        else
            B(kk+1,subevent) = B(kk,subevent) - 1;
            left_right = sign(rand-0.5);
            B(kk+1,subevent+left_right) = B(kk,subevent+left_right) + 1;
        end
    end
    %if sum(A(kk+1,:)) ~= Atot || sum(B(kk+1,:)) ~= Btot
    %    keyboard; % debug if molecules not conserved with pure diffusion
    %end
    time=time+tau;
    % plotting the process in real-time (downsampled for speed)
    if mod(kk,10000) == 0
        str = ['time = ',num2str(round(time)),' seconds'];
        figure(1);
        subplot(1,2,1), bar(A(kk,:));
        title('A molecules');
        xlabel('compartment number');
        ylim([0 500]);
        set(gca,'FontSize',20);
        grid on;
        subplot(1,2,2), bar(B(kk,:),'FaceColor',[0.9290 0.6940 0.1250]);
        ylim([0 120]);
        text(1,115,str,'FontSize',20);
        title('B molecules');
        xlabel('compartment number');
        set(gca,'FontSize',20);
        grid on;
    end
end
%% Plotting: Initial condition vs final state
figure;
subplot(1,2,1), bar(A(1,:));
title('Initial Condition');
ylabel('A molecules');
xlabel('x');
ylim([0 500]);
set(gca,'FontSize',20);
grid on;
subplot(1,2,2), bar(A(kk,:));
ylim([0 500]);
title('Final time');
xlabel('x');
set(gca,'FontSize',20);
grid on;

figure;
subplot(1,2,1), bar(B(1,:),'FaceColor',[0.9290 0.6940 0.1250]);
ylabel('B molecules');
xlabel('x');
title('Initial Condition');
ylim([0 120]);
set(gca,'FontSize',20);
grid on;
subplot(1,2,2), bar(B(kk,:),'FaceColor',[0.9290 0.6940 0.1250]);
ylim([0 120]);
xlabel('x');
title('Final time');
set(gca,'FontSize',20);
grid on;

%% Plotting: Final state plots

% str = ['time = ',num2str(round(time)),' seconds'];
% 
% figure;
% bar(A(kk,:));
% ylabel('A molecules');
% xlabel('x');
% text(10,480,str,'FontSize',20);
% ylim([0 500]);
% set(gca,'FontSize',20);
% grid on;
% 
% figure;
% bar(B(kk,:),'FaceColor',[0.9290 0.6940 0.1250]);
% ylim([0 120]);
% text(10,115,str,'FontSize',20);
% ylabel('B molecules');
% xlabel('x');
% set(gca,'FontSize',20);
% grid on;

%% Plotting: Space-time plots

figure;
imagesc(A(1:kk,:)), xlabel('compartment number'), ylabel('time (seconds)'),...
    set(gca,'YDir','normal');
title('A molecules');
clim([100 500]);
set(gca,'FontSize',20);
set(gca,'linewidth',2);
colorbar;
colormap jet;
yticks([0 floor(kk/4) floor(kk/2) floor(3*kk/4) kk]);
yticklabels({num2str(0), num2str(T/4), num2str(T/2), num2str(3*T/4), num2str(T)});

% figure;
% imagesc(B(1:kk,:)), xlabel('x'), ylabel('time (seconds)'),...
%     set(gca,'YDir','normal');
% title('B molecules');
% clim([50 80]);
% set(gca,'FontSize',20);
% set(gca,'linewidth',2);
% colorbar;
% colormap jet;
% yticks([0 floor(kk/4) floor(kk/2) floor(3*kk/4) kk]);
% yticklabels({num2str(0), num2str(T/4), num2str(T/2), num2str(3*T/4), num2str(T)});



