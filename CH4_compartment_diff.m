function CH4_compartment_diff

D=0.0001;         % 10^{-4} mm^2 sec^{-1}
h=0.025;

numberofrealisations=12;
Xinitial=0.4;

for i=1:numberofrealisations

    if (i<=numberofrealisations/2)
        X=Xinitial-h/2;
    else
        X=Xinitial+h/2;
    end
    time=0;
    k=1;
    Xplot(i,k)=X;
    timeplot(i,k)=time;

    while (k<320)
        rr=rand(2,1);
        a0=2*D/(h*h);
        time=time+(1/a0)*log(1/rr(1));
        if (rr(2)*a0<D/(h*h))
            X=X-h;
        else
            X=X+h;
        end
        if (X<0)
            X=h/2;
        end
        if (X>1)
            X=1-h/2;
        end
        k=k+1;
        Xplot(i,k)=X;
        timeplot(i,k)=time;
    end
end

timeplot=timeplot/60;

figure;
for i = 1:numberofrealisations
    stairs(Xplot(i,:),timeplot(i,:),'LineWidth',1.5);
    hold on;
end

xlabel('\eta_x','interpreter','tex');
ylabel('time [min]','interpreter','tex');
axis([-0.2 1.2 0 15]);
set(gca,'Fontsize',20);
grid on;

