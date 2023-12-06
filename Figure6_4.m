function Figure6_4

results=load('data_Figure6_4_a.dat');
resultsbeta0275=load('data_Figure6_4_b.dat');
Ms=0.25+(besseli(0,8)+besseli(2,8))/besseli(1,8);
set(groot,'defaultAxesTickLabelInterpreter','tex');

figure;
hh1=semilogx(1./[1 50],[Ms Ms],'m','Linewidth',4');
hold on;
hh3=semilogx(1./resultsbeta0275(:,1),resultsbeta0275(:,2),'ok','Markersize',4,'Linewidth',4);
semilogx([0.0275 0.0275],[2 3.3],'b--','Linewidth',4);
hh2=semilogx(1./results(:,1),results(:,2),'or','Markersize',4,'Linewidth',4);
hl=legend([hh1,hh2,hh3],'M_s (true mean)','\beta=0 (standard)','\beta=0.275 (rescaled)');
set(hl,'interpreter','tex','Fontsize',20,'location','northeast');
xlabel('h (log scale)','interpreter','tex');
ylabel('M(h)','interpreter','tex');
%set(gca,'XTick',[0.028 0.1 0.3 1]);
%set(gca,'XTickLabel',{'h_{\rm crit}' '10^{-1}' '3{}\times{}10^{-1}' '10^0'});
axis([0 1 2 3.2]);
set(gca,'Fontsize',20);
grid on;
