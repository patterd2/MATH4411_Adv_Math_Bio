k1 = 0.00025;
k2 = 0.18;
k3 = 37.5;

params = [1750 2100 2200 2450];
str = strings(length(params));

x = 50:1:450;
for i = 1:length(params)
    k4 = params(i);
    plot(x,-k1*x.^3 + k2*x.^2 -k3*x + k4,'LineWidth',3)
    hold on;
    xlabel('a');
    ylabel('f(a)');
    set(gca,'FontSize',20);
    str(i) = ['k_4 = ', num2str(k4)];
end
plot(x,zeros(1,length(x)),'-.k','LineWidth',2)
legend(str(1),str(2),str(3),str(4));
grid on;
axis tight;