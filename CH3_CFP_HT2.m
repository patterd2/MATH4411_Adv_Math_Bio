function CH3_CFP_HT2

load data_Figure3_4b_1.dat;
ys=data_Figure3_4b_1(:,1);
taus=data_Figure3_4b_1(:,2);

load data_Figure3_4b_2.dat;
y=data_Figure3_4b_2(:,1);
tau=data_Figure3_4b_2(:,2);

figure;
bar(ys,taus);
hold on;
axis([- 0.5 19 0 400]);
set(gca,'XTick',[0 2 4 6 8 10 12 14 16 18 20]);
set(gca,'YTick',[0 100 200 300 400 500 600 700]);
xlabel('y');
ylabel('E[\tau (y)] [sec]','interpreter','tex');
plot(y,tau,'r','Linewidth',3);
legend('stochastic simulation','integral formula');
set(gca,'Fontsize',20);
grid on;
