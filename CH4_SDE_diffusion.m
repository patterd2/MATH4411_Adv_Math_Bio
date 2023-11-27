function CH4_SDE_diffusion

finaltime=10*60;  % 10 minutes
D=0.0001;         % 10^{-4} mm^2 sec^{-1}

x=-1:0.05:1;
y=-1:0.05:1;
for i=1:41
    for j=1:41
        gauss(i,j)=(1/(4*D*pi*finaltime))*exp(-(x(i).*x(i)+y(j).*y(j))/(4*D*finaltime));
    end
end

figure(1);
pcolor(x,y,gauss');
title('P(\eta_x, \eta_y, t=10)');
xlabel('\eta_x','interpreter','tex');
ylabel('\eta_y','interpreter','tex');
set(gca,'XTick',[-0.6 -0.4 -0.2 0 0.2 0.4 0.6]);
set(gca,'YTick',[-0.6 -0.4 -0.2 0 0.2 0.4 0.6]);
axis([-0.7 0.7 -0.7 0.7]);
shading flat;
view(0,90);
box on;
colormap jet;
colorbar;
shading interp;
set(gca,'Fontsize',20);

dt=10;
n=finaltime/dt;
numberofrealizations=1000000;
hist=zeros(40,40);

for i=1:numberofrealizations
    dX = sqrt(2*D*dt)*randn(1,n);   
    X = cumsum(dX);             
    dY = sqrt(2*D*dt)*randn(1,n);   
    Y = cumsum(dY);
    if (((X(n)>-1)&&(X(n)<1))&&((Y(n)>-1)&&(Y(n)<1)))
            hist(round(20*X(n)+20.5),round(20*Y(n)+20.5))=hist(round(20*X(n)+20.5),round(20*Y(n)+20.5))+1;
    end
end

hist=hist/numberofrealizations/0.05/0.05;

xs=-1.:0.05:0.95;
ys=-1.:0.05:0.95;

figure;
pcolor(xs,ys,hist');
title('P(\eta_x, \eta_y, t=10) [SSA approx.]');
xlabel('\eta_x','interpreter','tex');
ylabel('\eta_y','interpreter','tex');
set(gca,'XTick',[-0.6 -0.4 -0.2 0 0.2 0.4 0.6]);
set(gca,'YTick',[-0.6 -0.4 -0.2 0 0.2 0.4 0.6]);
axis([-0.7 0.7 -0.7 0.7]);
shading faceted;
view(0,90);
box on;
colormap jet;
colorbar;
set(gca,'Fontsize',20);


