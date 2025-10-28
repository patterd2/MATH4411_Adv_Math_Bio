function CH4_SDE_reflect

dt=0.1;
finaltime=15*60;  % 15 minutes
D=0.0001;         % 10^{-4} mm^2 sec^{-1}

n=finaltime/dt+1;
numberofrealisations=10;
X=zeros(numberofrealisations,n);

for i=1:numberofrealisations
    X(i,1)=0.5;
    for j=2:n
        X(i,j)=X(i,j-1)+sqrt(2*D*dt)*randn(1,1);
        if (X(i,j)<0)
            X(i,j)=-X(i,j);
        end
        if (X(i,j)>1)
            X(i,j)=2-X(i,j);
        end
    end
end

time=(0:dt:finaltime)/60;

figure;
for i = 1:numberofrealisations
    plot(X(i,:),time);
    hold on;
end
xlabel('\eta_x','interpreter','tex');
ylabel('time [min]','interpreter','tex');
xlim([-0.2 1.2]);
set(gca,'Fontsize',20);
grid on;
box on;
