function CH1_naive_SSA_degradation

k=0.1; % reaction rate parameter
Ainitial=20; % initial number of A molecules, n0 in the notes
dt=0.005; % discretisation parameter/step size -> play with this, test limits!
numberofrealisations=50; % number of sample paths to run

times = 0:dt:30;
Adata = zeros(length(times),numberofrealisations);
%% Naive SSA of the degradation process
for i=1:numberofrealisations

    A=Ainitial;
    time=0;
    j=1;
    Aplot(1,i)=A;
    timeplot(1,i)=0;
    Adata(1,i) = A;
    m = 0;
    while (j<Ainitial+1)
        rr=rand(1,1);
        time=time+dt;
        m = m+1;
        if (rr<k*A*dt)
            A=A-1;
            j=j+1;
            Aplot(j,i)=A;
            timeplot(j,i)=time;
        end
        if m < length(times)+1
            Adata(m,i)=A;
        end
    end
    Aplot(j+1,i)=A;
    if (time<30)
        timeplot(j+1,i)=30;
    else
        timeplot(j+1,i)=time;
    end
end

tdet=zeros(301,1);
Adet=zeros(301,1);
tdet=[0:0.1:30]; % this is just for the deterministic comparison
Adet(:)=Ainitial*exp(-k*tdet(:));

%% Plot some sample paths for comparison to the mean (deterministic model)
figure; 
for i=1:min(2,numberofrealisations)
    h=stairs(timeplot(:,i),Aplot(:,i));
    set(h,'Linewidth',2);
    hold on;
end
plot(tdet,Adet,'k--','Linewidth',4);
xlabel('time [sec]');
ylabel('Number of Molecules');
if numberofrealisations > 2
    for i = 3:numberofrealisations
        h=stairs(timeplot(:,i),Aplot(:,i));
        set(h,'Linewidth',2);
        hold on;
    end
end
plot(tdet,Adet,'k--','Linewidth',4);
hh=legend('first realization','second realization','deterministic mean');
axis([0 30 0 (Ainitial+1)]);
set(gca,'linewidth',1.5);
set(gca,'FontSize',20);
grid on;
%% Plot the stochastic mean (estimated from multiple simulations) against the deterministic mean
%keyboard; % this is a useful command for debugging!
figure;
plot(tdet,Adet,'k--','Linewidth',4);
hold on;
plot(times,mean(Adata(:,:)'),'Linewidth',4);
xlabel('time [sec]');
ylabel('Number of Molecules');
str = ['stochastic mean (' num2str(numberofrealisations) ' realisations)'];
hh=legend('deterministic mean',str);
axis([0 30 0 (Ainitial+1)]);
set(gca,'linewidth',1.5);
set(gca,'FontSize',20);
grid on;

