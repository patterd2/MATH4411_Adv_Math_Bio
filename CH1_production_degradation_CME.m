function CH1_production_degradation_CME

k1=0.1; % production rate
k2=1; % degradation rate

numberofrealisations=100;
T = 100; % length of each simulation

M = 30; % keep track of A values between 0 and M
p=zeros(M,1); % need to adjust size when parameters are adjusted

for i = 1:numberofrealisations
    time=0;
    A = 15; %ceil(k2/k1); % initialise at the expected mean to save time
    while (time < T)
        rr=rand(2,1);
        a0=k1*A+k2;
        time=time+(1/a0)*log(1/rr(1));
        if (rr(2)*a0<k1*A)
            A=A-1; % degradation
        else
            A=A+1; % production
        end
    end
    if A < M
    p(A+1)=p(A+1)+1; % keep track of num molecules at simulation end
    end
end
p=p/numberofrealisations;

%% Compute the stationary distribution using the recursive formulation
% p(n+1) is the stationary probability that A=n for n=0,1,2,3,.....

phi(1)=1;
phi(2)=k2/k1*phi(1);
int=phi(1)+phi(2);
for n=1:24
    phi(n+2)=(k1*n*phi(n+1)+k2*phi(n+1)-k2*phi(n))/(k1*(n+1));
    int=int+phi(n+1);
end
phi=phi/int; % impose normalisation constraint (must be a true PMF!)

%% Plotting
figure;
bar(0:23,p(1:24));
hold on;
plot(0:25,phi,'r','Linewidth',3);
xlabel('number of molecules','Interpreter','none');
ylabel('probability','Interpreter','none');
legend('Estimated PMF from SSA','Stationary distribution','Interpreter','none');
axis tight;
xlim([0 25]);
set(gca,'linewidth',1.5);
set(gca,'FontSize',20);
grid on;
