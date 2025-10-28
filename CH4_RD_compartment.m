% second sim
function CH4_RD_compartment

h=0.01;
xdet=[h/2:h:1-h/2];
w=size(xdet);
N=w(2);

[t,y]=ode45(@density,[0 2000*60],[zeros(N,1) zeros(N,1)]);
z=size(y);
endline=z(1);
detdenA=y(endline,1:N)';
detdenB=y(endline,N+1:2*N)';

histogram=load('data_Figure6_3.dat');
ww=size(histogram);
boxes=ww(1);
mesh=[1/(2*boxes):1/boxes:1-1/(2*boxes)]';

figure;
h1=bar(mesh,histogram(:,1),'b');
hold on;
h2=plot([0 xdet 1],[detdenA(1); detdenA; detdenA(N)],'r','Linewidth',4);
box on;
xlabel('x','interpreter','tex');
ylabel('number of A molecules','interpreter','tex');
hh=legend([h1 h2],'Gillespie SSA','Mass action PDE');
axis([0 1 0 18]);
set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1]);
set(gca,'Fontsize',20);
grid on;

figure;
h1=bar(mesh,histogram(:,2),'b');
hold on;
h2=plot([0 xdet 1],[detdenB(1); detdenB; detdenB(N)],'r','Linewidth',4);
box on;
xlabel('x','interpreter','tex');
ylabel('number of B molecules','interpreter','tex');
hh=legend([h1 h2],'Gillespie SSA','Mass action PDE');
axis([0 1 0 120]);
set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1]);
set(gca,'Fontsize',20);
grid on;


function dy=density(t,y)
h=0.01;
N=100;
D=0.000001;
Diff=D/(h*h);
hbox=0.025;
k1=0.03/(1000*1000*1000)/(hbox*hbox*hbox);
k2=0.3/(1000*1000*1000)/(hbox*hbox*hbox);
k3=0.0001;
k4=0.0001;
k5=0.0000001*1000*1000*1000*(hbox*hbox*hbox);
k6=0.000001*1000*1000*1000*(hbox*hbox*hbox);

Dmatrix=-2*diag(ones(N,1),0)+diag(ones(N-1,1),1)+diag(ones(N-1,1),-1);
Dmatrix(1,1)=-1;
Dmatrix(N,N)=-1;
Dmatrix=Diff*Dmatrix;

dy=zeros(2*N,1);
dy(1:N)=Dmatrix*y(1:N)-2*k1*y(1:N).*y(1:N)-k2*y(1:N).*y(N+1:2*N)+k5*ones(N,1)-k3*y(1:N);
dy(N+1:2*N)=Dmatrix*y(N+1:2*N)-k2*y(1:N).*y(N+1:2*N)-k4*y(N+1:2*N)+k6*[zeros(3*N/5,1); ones(2*N/5,1)];
