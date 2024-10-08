%% Plot the deterministic direction field for part (a)
x = -0.25:0.01:1.25;
f = -x.*(x-0.5).*(x-1);
figure;
plot(x,f,'LineWidth',3);
xlabel('x');
ylabel('f(x)');
hold on;
tempX = -0.2:0.1:0;
plot(tempX,zeros(1,length(tempX)),'>m','LineWidth',2); 
tempX = 0:0.1:0.5;
plot(tempX,zeros(1,length(tempX)),'<m','LineWidth',2); 
tempX = 0.5:0.1:1;
plot(tempX,zeros(1,length(tempX)),'>m','LineWidth',2); 
tempX = 1:0.1:1.2;
plot(tempX,zeros(1,length(tempX)),'<m','LineWidth',2); 
scatter(0,0,100,'red','filled');
scatter(0.5,0,100,'black','filled');
scatter(1,0,100,'red','filled');
ylim([-0.06 0.06]);
xlim([-0.25 1.25]);

%% Plot the stationary distribution for part (b)
x = -0.5:0.01:1.5;
ps = exp( -25*x.^2 + 50*x.^3 -25*x.^4 );
ps = ps./trapz(x,ps);
figure;
plot(x,ps,'LineWidth',4);
xlabel('x');
ylabel('p_s(x)','Interpreter','tex');
axis tight;
hold on;
scatter(0,0,100,'red','filled');
scatter(0.5,0,100,'black','filled');
scatter(1,0,100,'red','filled');

%% Plot the stationary distribution for part (c)
x = -1:0.01:1;
ps = exp(4*x.^4);
ps = ps./trapz(x,ps);
figure;
plot(x,ps,'b','LineWidth',4);
xlabel('x');
ylabel('p_{s,2}(x)','Interpreter','tex');
hold on;
x2 = -3:0.01:-1;
plot(x2,zeros(1,length(x2)),'b','LineWidth',4);
x2 = 1:0.01:3;
plot(x2,zeros(1,length(x2)),'b','LineWidth',4);
xline(1,'--b','LineWidth',2)
xline(-1,'--b','LineWidth',2)
axis tight;

%%
x = -1:0.01:1;
ps = exp(4*x.^4);
ps = ps./trapz(x,ps);
figure;
plot(x,ps,'b','LineWidth',4);
xlabel('x');
ylabel('p_{s,3}(x)','Interpreter','tex');
hold on;
x2 = -3:0.01:-1;
plot(x2,ps,'b','LineWidth',4);
x2 = 1:0.01:3;
plot(x2,ps,'b','LineWidth',4);
axis tight;

%% Problem Class 5, Q1 (b)
x = 0:0.01:1;
delta = 0.1;
beta = 0.75;
phi = 0.75;
omega = 0.2;
plot(x,delta-beta*x+phi*x.^2);
hold on;
G1 = (beta+sqrt(beta^2 - 4*delta*phi))/(2*phi);
G2 = (beta-sqrt(beta^2 - 4*delta*phi))/(2*phi);
scatter(G1,0,'filled');
scatter(G2,0,'filled');

T1 = (omega*(1-G1))/(delta+phi*G1*G1+omega);
T2 = (omega*(1-G2))/(delta+phi*G2*G2+omega);

J_10 = [0 delta+phi-beta; -omega -delta-phi-omega];

J_G1_T1 = [2*phi*G1*T1-beta*T1 delta+phi*G1^2-beta*G1; -2*phi*G1*T1-omega -delta-phi*G1^2-omega];
J_G2_T2 = [2*phi*G2*T2-beta*T2 delta+phi*G2^2-beta*G2; -2*phi*G2*T2-omega -delta-phi*G2^2-omega];