%using optimized cmaes code and modifying it for the current application
tic()
%initialization
N=25;
lamda=25;
stopfitness=0.1;
stopeval=10*(1e3*N^2);
%strategy parameters: selection
lambda=4+floor(3*log(N));
mu=lambda/2;
weights=log(mu+(1/2))-log(1:mu)';
mu=floor(mu);
weights=(weights/sum(weights));
mueff=(sum(weights)^2)/(sum(weights.^2));
%strategy parameters:Adaption
cv=2/((N+1.3)^2+mueff);
cc=(4+mueff/N)/(N+(4+2*mueff/N));
csig=(mueff+2)/(N+mueff+5);
cmu=2*(mueff-2+(1/mueff))/((N+2)^2+(2*(mueff/2)));
damp=1+2*max(0,sqrt((mueff-1)/(N+1))-1)+csig;
pc=zeros(N,1);psig=zeros(N,1);sig=0.5;
m=randn(N,1);
Q=zeros(1,1000000);
t=0:0.01:7;
 [x1,x2,x3,x4,x5,x6,u,v,x,G,U,En1,En2,En3,y]=coupled_doubleinv_exp(m);
 x1i=x1;
 x2i=x2;
 x3i=x3;
 x4i=x4;
 x5i=x5;
 x6i=x6;
Q(1)=y;
r=zeros(N,lamda);
I=eye(N);
B=eye(N);
D=eye(N);
C=B*D*(B*D)';
chiN=N^0.5*(1-1/(4*N)+1/(21*N^2));
eigenval=0;
counteval=0;
Em=zeros(lamda,1);
epsi=zeros(N,lamda);
for G1=1:1000000
for i=1:lamda
    epsi(:,i)=randn(N,1);
    r(:,i)=m+sig*(B*D*epsi(:,i));
    [x1,x2,x3,x4,x5,x6,u,v,x,G,U,En1,En2,En3,y]=coupled_doubleinv_exp(r(:,i));
    Em(i)=y;
    counteval=counteval+1;
end
[E2,I]=sort(Em);
m=r(:,I(1:mu))*weights;
m1=epsi(:,I(1:mu))*weights;
y1=E2(1)
Q(G1)=E2(1);
%cumumation update
psig=(1-csig)*psig+(sqrt(csig*(2-csig)*mueff))*(B*m1);
hsig=norm(psig)/sqrt(1-(1-csig)^(2*counteval/lamda))/chiN;
pc=(1-cc)*pc+hsig*(sqrt(cc*((2-cc))*mueff))*(B*D*m1);
C=(1-cv-cmu)*C+cv*(pc*pc'+(1-hsig)*cc*(2-cc)*C)+cmu*(B*D*epsi(:,I(1:mu)))*diag(weights)*(B*D*epsi(:,I(1:mu)))';
%step_size_adaption 
sig=sig*exp((csig/damp)*(((norm(psig))/(chiN))-1));
%update B and D
if counteval-eigenval>lamda/(cv+cmu)/N/10
eigenval=counteval;
C=triu(C)+triu(C,1)';
[B,D]=eig(C);
D=diag(sqrt(diag(D)));
end
if (E2(1)<=stopfitness)
    break;
end
if E2(1)==E2(ceil(0.7*lamda))
    sig=sig*exp(0.2+csig/damp);
end
plot(t,U)
pause(0.01);
end
 [x1,x2,x3,x4,x5,x6,u,v,x,G,U,En1,En2,En3,y]=coupled_doubleinv_exp(m);
toc()

figure;plot(t,U);title('control action');xlabel('time (in seconds)','interpreter','latex','FontWeight','bold');ylabel('$u$ (in Newtons)','interpreter','latex','FontWeight','bold');grid on;
figure;semilogx(Q(1:G1));ylabel('$f_1$','interpreter','latex','FontWeight','bold');xlabel('Generation (logarithmic scale)','FontWeight','bold');grid on;axis('tight');
figure;subplot(1,2,1);plot(t,x1i);title('before');ylabel('$y$ (in meters)','interpreter','latex','FontWeight','bold');xlabel('time (in seconds)','interpreter','latex','FontWeight','bold');axis('tight');grid on;subplot(1,2,2);plot(t,x1);title('after');ylabel('$y$ (in meters)','interpreter','latex','FontWeight','bold');xlabel('time (in seconds)','interpreter','latex','FontWeight','bold');grid on;axis('tight');
figure;subplot(1,2,1);plot(t,cos(x2i));title('before');xlabel('time (in seconds)','interpreter','latex','FontWeight','bold');ylabel('cos($\theta_1$)','interpreter','latex','FontWeight','bold');grid on;axis('tight');subplot(1,2,2);plot(t,cos(x2));title('after');xlabel('time (in seconds)','interpreter','latex','FontWeight','bold');ylabel('cos($\theta_1$)','interpreter','latex','FontWeight','bold');axis('tight');grid on;
figure;axis('tight');subplot(1,2,1);plot(t,cos(x3i));title('before');xlabel('time (in seconds)','interpreter','latex','FontWeight','bold');ylabel('cos($\theta_2$)','interpreter','latex','FontWeight','bold');axis('tight');grid on;subplot(1,2,2);plot(t,cos(x3));title('after');xlabel('time (in seconds)','interpreter','latex','FontWeight','bold');ylabel('cos($\theta_2$)','interpreter','latex','FontWeight','bold');grid on;axis('tight');
figure;subplot(1,2,1);plot(t,x4i);title('before');xlabel('time (in seconds)','interpreter','latex','FontWeight','bold');ylabel('$\dot{y}$ (in meters/second)','interpreter','latex','FontWeight','bold');grid on;axis('tight');subplot(1,2,2);plot(t,x4);title('after');xlabel('time (in seconds)','interpreter','latex','FontWeight','bold');ylabel('$\dot{y}$ (in meters/second)','interpreter','latex','FontWeight','bold');grid on;axis('tight');
figure;subplot(1,2,1);plot(t,x5i);title('before');grid on;xlabel('time (in seconds)','interpreter','latex','FontWeight','bold');ylabel('$\dot{\theta_1}$ (in radians/second)','interpreter','latex','FontWeight','bold');grid on;axis('tight');subplot(1,2,2);plot(t,x5);title('after');xlabel('time (in seconds)','interpreter','latex','FontWeight','bold');ylabel('$\dot{\theta_1}$ (in radians/second)','interpreter','latex','FontWeight','bold');axis('tight');grid on;
figure;subplot(1,2,1);plot(t,x6i);title('before');grid on;xlabel('time (in seconds)','interpreter','latex','FontWeight','bold');ylabel('$\dot{\theta_2}$ (in radians/second)','interpreter','latex','FontWeight','bold');axis('tight');subplot(1,2,2);plot(t,x6);title('after');xlabel('time (in seconds)','interpreter','latex','FontWeight','bold');ylabel('$\dot{\theta_2}$ (in radians/second)','interpreter','latex','FontWeight','bold');grid on;axis('tight');
figure;plot(t,u);xlabel('time (in seconds)','interpreter','latex','FontWeight','bold');ylabel('$x$ (in meters)','interpreter','latex','FontWeight','bold');axis('tight');grid on;
figure;plot(t,v);xlabel('time (in seconds)','interpreter','latex','FontWeight','bold');ylabel('$\dot{x}$ (in meters/second)','interpreter','latex','FontWeight','bold');grid on