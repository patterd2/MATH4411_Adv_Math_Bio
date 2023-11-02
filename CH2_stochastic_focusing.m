function CH2_stochastic_focusing
%% set parameters
k1=100;
k2=1000;
k3=0.01;
k4=9900;
k5=1000;
k6=100;

numberofrealisations = 3;

A0 = 10;
B0 = 100;
C0 = 0;
%% SSA
for i=1:numberofrealisations
    A=A0;
    B=B0;
    C=C0;
    time=0;
    for kk=1:2100
        while (time<kk*1)
            if (time>600)
                k5=500;
            end
            rr=rand(2,1);
            a0=k1+k2*C+k3*B+k4*A*C+k5+k6*A;
            tau=(1/a0)*log(1/rr(1));
            ss=k1;
            if (ss>rr(2)*a0)
                C=C+1;
            else
                ss=ss+k2*C;
                if (ss>rr(2)*a0)
                    C=C-1;
                    B=B+1;
                else
                    ss=ss+k3*B;
                    if (ss>rr(2)*a0)
                        B=B-1;
                    else
                        ss=ss+k4*A*C;
                        if (ss>rr(2)*a0)
                            C=C-1;
                        else
                            ss=ss+k5;
                            if (ss>rr(2)*a0)
                                A=A+1;
                            else
                                ss=ss+k6*A;
                                if (ss>rr(2)*a0)
                                    A=A-1;
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

[t,z] = ode15s(@myode,[0 2100],[10; 100; 0]);

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
axis([0 35 0 22]);
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
axis([0 35 50 400]);
set(gca,'Fontsize',20);
grid on;

%%
figure(3);
plot(t(1:3:end),z(1:3:end,3),'-.k','Linewidth',2);
hold on;
for j = 1:1
    [~, t2] = max(tt(:,j));
    scatter(tt(1:t2,j),CC(1:t2,j),'Linewidth',4);
end
plot(t(1:3:end),z(1:3:end,3),'-.k','Linewidth',2);
xlabel('time [min]');
ylabel('number of C molecules');
legend('solution of ODEs');
axis([0 35 0 1.2]);
set(gca,'Fontsize',20);
grid on;

averA=10;
detB=(k2/k3)*k1/(k2+k4*averA);
n=0:50;
stoB=sum(((k1*k2/k3)./(k2+k4*n)).*(((averA).^n)./factorial(n))*exp(-averA));
disp('time < 10 min');
display(['deterministic B = ',num2str(detB)]);
display(['stochastic B = ',num2str(stoB)]);

averA=5;
detB=(k2/k3)*k1/(k2+k4*averA);
stoB=sum(((k1*k2/k3)./(k2+k4*n)).*(((averA).^n)./factorial(n))*exp(-averA));
disp('time > 10 min');
disp(['deterministic B = ',num2str(detB)]);
disp(['stochastic B = ',num2str(stoB)]);


%% Define RHS of ODE
function dydt = myode(t,z)
k1=100;
k2=1000;
k3=0.01;
k4=k2*99/10;
k5=1000;
k6=100;
if (t>600)
    k5=500;
end
dydt =[k5-k6*z(1); k2*z(3)-k3*z(2); k1-k2*z(3)-k4*z(1)*z(3)];