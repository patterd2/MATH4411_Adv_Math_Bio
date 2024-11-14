%function PC3_Q2_stochastic_focusing
tic;
%% set parameters
k1 = 300;
k2 = 20000;
k3 = 1;
k4 = 500;
k5 = 100;
k6 = 200000;

% Check if we are in the right parameter regime
%display((k5^2 + k1*k5)/(k2*k4 + k5*k6));

numberofrealisations = 250;

% A0 = 10;
% B0 = 100;
% C0 = 0;

A0 = ceil(k4/k5);
C0 = 0;
B0 = ceil(k2*A0*C0/k3);

%% SSA
for i=1:numberofrealisations
    k4 = 500;
    A=A0;
    B=B0;
    C=C0;
    time=0;
    for kk=1:2100
        while (time<kk*1)
            if (time>600)
                k4=1000;
            end
            rr=rand(2,1);
            a0=k1+k2*A*C+k3*B+k4+k5*A+k6*C;
            tau=(1/a0)*log(1/rr(1));
            ss=k1;
            if (ss>rr(2)*a0)
                C=C+1;
            else
                ss=ss+k2*A*C;
                if (ss>rr(2)*a0)
                    C=C-1;
                    B=B+1;
                else
                    ss=ss+k3*B;
                    if (ss>rr(2)*a0)
                        B=B-1;
                    else
                        ss=ss+k4;
                        if (ss>rr(2)*a0)
                            A=A+1;
                        else
                            ss=ss+k5*A;
                            if (ss>rr(2)*a0)
                                A=A-1;
                            else
                                ss=ss+k6*C;
                                if (ss>rr(2)*a0)
                                    C=C-1;
                                end
                            end
                        end
                    end
                end
            end
            time=time+tau;
        end
        AA(kk,i)=A;
        BB(kk,i)=B;
        CC(kk,i)=C;
        tt(kk,i)=time/60;
    end
end

[t,z] = ode15s(@myode,[0 2200],[A0; B0; C0]);

t=t/60;
%% Plotting
figure(1);
plot(t,z(:,1),'k','Linewidth',4);
hold on;
for j = 1:numberofrealisations
    [~, t2] = max(tt(:,j));
    stairs(tt(1:t2,j),AA(1:t2,j),'Linewidth',2);
end
plot(t,z(:,1),'k','Linewidth',4);
xlabel('time [min]');
ylabel('number of A molecules');
legend('solution of ODEs');
axis tight;
set(gca,'Fontsize',20);
grid on;

%%
figure(2);
plot(t,z(:,2),'k','Linewidth',4);
hold on;
for j = 1:numberofrealisations
    [~, t2] = max(tt(:,j));
    stairs(tt(1:t2,j),BB(1:t2,j),'Linewidth',2);
end
plot(t,z(:,2),'k','Linewidth',4);
xlabel('time [min]');
ylabel('number of B molecules');
legend('solution of ODEs');
axis tight;
set(gca,'Fontsize',20);
grid on;

%% Compute the approximations for the mean of B
sum = 0;
for n = 0:50
    sum = sum + (1/factorial(n))*(k2*n/k3)*(k1/(k6+k2*n))*exp(-k4/k5)*(k4/k5)^n ;
end
B_bar = (k1*k2*k4)/(k3*(k5*k6 + k2*k4));
disp(['The monte-carlo approx. of the mean of B is: ',num2str(mean(BB(:,end)))]);
disp(['The stochastic approx. of the mean of B is: ',num2str(sum)]);
disp(['The mass action approx. of the mean of B is: ',num2str(B_bar)]);
%%
figure(3);
plot(t(1:3:end),z(1:3:end,3),'-.k','Linewidth',3);
hold on;
for j = 1:1
    [~, t2] = max(tt(:,j));
    scatter(tt(1:t2,j),CC(1:t2,j),'Linewidth',4);
end
plot(t(1:3:end),z(1:3:end,3),'-.k','Linewidth',3);
xlabel('time [min]');
ylabel('number of C molecules');
legend('solution of ODEs');
axis tight;
set(gca,'Fontsize',20);
grid on;

toc;
%% Define RHS of ODE
function dydt = myode(t,z)
k1 = 300;
k2 = 20000;
k3 = 1;
k4 = 500;
k5 = 100;
k6 = 200000;

if (t>600)
    k4 = 1000;
end
dydt =[k4-k5*z(1); k2*z(1)*z(3)-k3*z(2); k1-k6*z(3)-k2*z(1)*z(3)];
end