% this code does phase plane analysis for the SIR epidemic model
function SIR_pp
global alp bet s0

set(0,                           ...
    'defaultaxesfontsize', 20,   ...
    'defaultaxeslinewidth', 1.0, ...
    'defaultlinelinewidth', 1.2, ...
    'defaultpatchlinewidth', 0.7);

%set parameters
alp = 1;
bet = 1;

Imx = 1;% limit on y-axis plots

%nullcline

% now create some phase portraits
tstep = 0.01; %time step size
t_end = 50; % length of time interval
tspan = [0:tstep:t_end];

%initial data
s0_list = [1.5,2,2.5];
for j = 1:3

    s0 = [s0_list(j),0.001,0]; % initial data for integration

    [T,S] = ode23s(@deRHS,tspan,s0);
    figure(1)
    plot(S(:,1),S(:,2),'linewidth',2)
    xlabel('\alpha s/\beta','fontsize',20)
    ylabel('\alpha i/\beta','fontsize',20)
    axis([0 3 0 Imx])
    hold on
end

figure(1);
plot([bet/alp,bet/alp],[0,Imx],'--','linewidth',2);
annotation('textarrow',[0.42 .35],[0.68 0.68],'linewidth',2);
hold on;
[X,Y] = meshgrid(0.2:0.2:3,0.2:0.2:3);
du = -alp*X.*Y;
dv = alp*X.*Y-bet*Y;
un=du./sqrt(du.^2+dv.^2);
uv=dv./sqrt(du.^2+dv.^2);
quiver(X,Y,un,uv,0.5);
grid on;
hold off

%print('../../figs_c/chapt_1/sir_pp','-deps','-cmyk')
sinfbys0 = [0.01:0.01:1];
as0byb = log(sinfbys0)./(sinfbys0-1);


figure(2)
plot(as0byb,sinfbys0,'linewidth',2)

xlabel('R_0=\alpha s(0)/\beta','fontsize',20)
ylabel('s(\infty)/s(0)','fontsize',20)
grid on;


function s_prime=deRHS(t,s)  % right hand side for ode system
global alp bet

S = s(1);
I = s(2);
R = s(3);

dsdt = -alp*S*I;
didt = alp*S*I-bet*I;
drdt = bet*I;

s_prime = [dsdt didt drdt]';
