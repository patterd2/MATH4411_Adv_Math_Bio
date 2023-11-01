function CH2_stochastic_resonance

%% Set parameter values
k1=0.0004;
k2=50;
k3=10;
k4=25;

A0=10;
B0=10;

time=0;
timeSSA=0;
kk=0;
A = A0;
B = B0;
%% SSA
while (time<500)
    timeSSAprev=timeSSA;
    timefinSSA=500;
    kk=kk+1;
    while (timeSSA<timeSSAprev+timefinSSA)
        if (A>100)
            timefinSSA=10;
        end
        rr=rand(2,1);
        timeSSA=timeSSA+1;
        a0=k1*A*(A-1)*B+k2+k3*A+k4;
        tau=(1/a0)*log(1/rr(1));
        ss=k1*A*(A-1)*B;
        if (ss>rr(2)*a0)
            A=A+1;
            B=B-1;
        else
            ss=ss+k2;
            if (ss>rr(2)*a0)
                A=A+1;
            else
                ss=ss+k3*A;
                if (ss>rr(2)*a0)
                    A=A-1;
                else
                    ss=ss+k4;
                    if (ss>rr(2)*a0)
                        B=B+1;
                    end
                end
            end
        end
        time=time+tau;
    end
    AA(kk)=A;
    BB(kk)=B;
    tt(kk)=time;
end

%% Solve ODE (mass action model)
[t,z] = ode45(@(t,z) [k1*z(1)*z(1)*z(2)+k2-k3*z(1); -k1*z(1)*z(1)*z(2)+k4] ,[0 500],[A0; B0]);

%% Plotting
figure(1); % logged y-axis version of the plot
h=semilogy(tt/60,AA);
hold on;
set(h,'Color','b','Linewidth',2);
h=semilogy(t/60,z(:,1));
set(h,'Color','r','Linewidth',4);
xlabel('time [min]');
ylabel('number of A molecules');
legend('stochastic','deterministic');
axis tight;
set(gca,'YTick',[1 10 100 1000 10000]);
set(gca,'Fontsize',20);

%%
figure(2);
plot(tt/60,AA,'r--','Linewidth',4);
hold on;
plot(tt/60,BB,'b--','Linewidth',2);
xlabel('time [min]');
ylabel('number of molecules');
hh=legend('A(t)','B(t)');
set(hh,'interpreter');
axis tight;
grid on;
set(gca,'Fontsize',20);

%%
figure(3);
xval1=[2:0.15:24.8 25:1:100 110:10:1000 1100:100:6200];
xval2=[2:0.15:24.8 25:1:160];
ynul1=(k3*xval1-k2)./(k1*xval1.*xval1);
ynul2=k4./(k1*xval2.*xval2);
semilogx(AA,BB,'b','Linewidth',2);
hold on;
h=semilogx(z(1,1),z(1,2),'b','MarkerSize',6);
set(h,'Color','r','Linewidth',2);
xlabel('number of A molecules');
ylabel('number of B molecules');
h=semilogx(xval1,ynul1);
set(h,'Color','g','Linewidth',4);
h=semilogx(xval2,ynul2);
set(h,'Color','g','Linewidth',4);
h=semilogx([(k4+k2)/k3],[k4*k3*k3/(k1*(k4+k2)*(k4+k2))],'ro');
set(h,'Markersize',4,'Linewidth',6);
h=semilogx(z(:,1),z(:,2));
set(h,'Color','r','Linewidth',2);
legend('stochastic','deterministic');
set(gca,'XTick',[1 10 100 1000 10000]);
set(gca,'Fontsize',20);
axis tight;
ylim([0 1500]);
xlim([0 1000]);

%%
figure(4);
[t,z] = ode45(@(t,z) [k1*z(1)*z(1)*z(2)+k2-k3*z(1); -k1*z(1)*z(1)*z(2)+k4] ,[0 500],[A0; B0]);
xval1=1:500;
xval2=1:1500;
ynul1=(k3*xval1-k2)./(k1*xval1.*xval1);
ynul2=k4./(k1*xval2.*xval2);
plot(xval1,ynul1,'g','Linewidth',4);
hold on;
plot(xval2,ynul2,'g','Linewidth',4);
ylim([0 1500]);
xlim([0 200]);
xlabel('number of A molecules');
ylabel('number of B molecules');
plot(z(:,1),z(:,2),'r','Linewidth',2);
scatter(z(1,1),z(1,2),'b','Linewidth',6);
set(gca,'Fontsize',20);
grid on;

