function CH3_example_2

dt=0.1;
finaltime=10*60;  % 10 minutes
D=0.0001;         % 10^{-4} mm^2 sec^{-1}
numberofrealisations=5;

n=finaltime/dt;
X=zeros(numberofrealisations,n);
Y=zeros(numberofrealisations,n);
Z=zeros(numberofrealisations,n);

for i=1:numberofrealisations
    dX = sqrt(2*D*dt)*randn(1,n);
    X(i,:)=cumsum(dX(1,:));
    dY = sqrt(2*D*dt)*randn(1,n);
    Y(i,:) = cumsum(dY(1,:));
    dZ = sqrt(2*D*dt)*randn(1,n);
    Z(i,:) = cumsum(dZ(1,:));
end

% need to downsample for readability
Xplot(:,:)=X(:,5:5:n);
Yplot(:,:)=Y(:,5:5:n);
Zplot(:,:)=Z(:,5:5:n);


%% Plotting
figure(1);
for j = 1:numberofrealisations
    plot([0,Xplot(j,:)],[0,Yplot(j,:)]);
    hold on;
    plot([X(j,n)],[Y(j,n)],'ko','Linewidth',5);
    plot(0,0,'ro','Linewidth',5);
end
xlabel('X [mm]');
ylabel('Y [mm]');
axis tight
set(gca,'Fontsize',20);
grid on;

figure(2);
for j = 1:min(1,numberofrealisations)
    plot3([0,Xplot(j,:)],[0,Yplot(j,:)],[0,Zplot(j,:)])
    hold on;
    plot3([X(j,n)],[Y(j,n)],[Z(j,n)],'ko','Linewidth',5);
    plot3(0,0,0,'ro','Linewidth',5);
end
xlabel('X [mm]');
ylabel('Y [mm]');
zlabel('Z [mm]');
axis tight
set(gca,'Fontsize',20);
grid on;

