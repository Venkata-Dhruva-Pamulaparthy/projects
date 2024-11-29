function[x1,x2,x3,x4,x5,x6,u,v,x,G,U,En1,En2,En3,y]=coupled_doubleinv_exp(m)
x1=zeros(1,701);
x2=zeros(1,701);
x3=zeros(1,701);
x4=zeros(1,701);
x5=zeros(1,701);
x6=zeros(1,701);
dt=0.01;
t=0:dt:7;
x1(1)=0;
x2(1)=pi;
x3(1)=pi;
x4(1)=0;
x5(1)=0;
x6(1)=0;
p1=m';
alpha=5;
beta=alpha/4;
gama=alpha/3;
u=zeros(length(t),1);
v=zeros(length(t),1);
x=zeros(length(t),1);
G=zeros(25,length(t));
U=zeros(1,length(t));
x0=1;
u0=1;
v0=0;
u(1)=u0;
v(1)=v0;
x(1)=x0;
%parameters
m=0.3;M=1;
l1=1;
l2=2;
g=9.8;
mu0=5;
tau=1;
c=[1.0000;0.7642;0.5840;0.4463;0.3411;0.2607;0.1992;0.1522;0.1163;0.0889;0.0679;0.0519;0.0397;0.0303;0.0232;0.0177;0.0135;0.0103;0.0079;0.0060;0.0046;0.0035;0.0027;0.0021; 0.0016];

h=25*1.5*[1/(c(1)^2);1/(c(2)^2);1/(c(3)^2);1/(c(4)^2);1/(c(5)^2);1/(c(6)^2);1/(c(7)^2);1/(c(8)^2);1/(c(9)^2);1/(c(10)^2);1/(c(11)^2);1/(c(12)^2);1/(c(13)^2);1/(c(14)^2);1/(c(15)^2);1/(c(16)^2);1/(c(17)^2);1/(c(18)^2);1/(c(19)^2);1/(c(20)^2);1/(c(21)^2);1/(c(22)^2);1/(c(23)^2);1/(c(24)^2);1/(c(25)^2)];

for i=1:length(t)-1
%parameter matrix
M11=M+2*m;
M12=-m*l1*cos(x2(i));
M13=-m*l2*cos(x3(i));
M21=-m*l1*cos(x2(i));
M22=m*(l1^2);
M23=0;
M31=-m*l2*cos(x3(i));
M32=0;
M33=m*(l2^2);
%nonlinearities
N11=(m*11*sin(x2(i))*(x5(i)^2))+(m*l2*sin(x3(i))*(x6(i)^2));
N21=-m*g*l1*sin(x2(i));
N31=-m*g*l2*sin(x3(i));

%inverse
detM=M11*(M22*M33-M32*M23)-M12*(M21*M33-M31*M23)+M13*(M21*M32-M31*M22);
MI11=(1/detM)*(M22*M33-M32*M23);
MI12=-(1/detM)*(M12*M33-M32*M13);
MI13=(1/detM)*(M12*M23-M22*M13);
MI21=-(1/detM)*(M21*M33-M31*M23);
MI22=(1/detM)*(M11*M33-M31*M13);
MI23=-(1/detM)*(M11*M23-M21*M13);
MI31=(1/detM)*(M21*M32-M31*M22);
MI32=-(1/detM)*(M11*M32-M31*M12);
MI33=(1/detM)*(M11*M22-M21*M12);

%state equations (Full dynamics of a parallel double inverted pendulum with a spring-like controller)
G(:,i)=[exp(-h(1)*((x(i)-c(1))^2));exp(-h(2)*((x(i)-c(2))^2));exp(-h(3)*((x(i)-c(3))^2));exp(-h(4)*((x(i)-c(4))^2));exp(-h(5)*((x(i)-c(5))^2));exp(-h(6)*((x(i)-c(6))^2));exp(-h(7)*((x(i)-c(7))^2));exp(-h(8)*((x(i)-c(8))^2));exp(-h(9)*((x(i)-c(9))^2));exp(-h(10)*((x(i)-c(10))^2));exp(-h(11)*((x(i)-c(11))^2));exp(-h(12)*((x(i)-c(12))^2));exp(-h(13)*((x(i)-c(13))^2));exp(-h(14)*((x(i)-c(14))^2));exp(-h(15)*((x(i)-c(15))^2));exp(-h(16)*((x(i)-c(16))^2));exp(-h(17)*((x(i)-c(17))^2));exp(-h(18)*((x(i)-c(18))^2));exp(-h(19)*((x(i)-c(19))^2));exp(-h(20)*((x(i)-c(20))^2));exp(-h(21)*((x(i)-c(21))^2));exp(-h(22)*((x(i)-c(22))^2));exp(-h(23)*((x(i)-c(23))^2));exp(-h(24)*((x(i)-c(24))^2));exp(-h(25)*((x(i)-c(25))^2))];
v(i+1)=v(i)+dt*((((mu0^2)*(0-u(i)))-2*tau*mu0*v(i))+(((p1*G(:,i))/(sum(G(:,i))))*(x(i)*(0-u0))));
u(i+1)=u(i)+dt*v(i);
U(i)=((((mu0^2)*(0-u(i)))-2*tau*mu0*v(i))+(((p1*G(:,i))/(sum(G(:,i))))*(x(i)*(0-u0))));
x(i+1)=x(i)-dt*(1/5)*(gama*x(i));
x1(i+1)=x1(i)+dt*(x4(i));
x2(i+1)=x2(i)+dt*(x5(i));
x3(i+1)=x3(i)+dt*(x6(i));
x4(i+1)=x4(i)+dt*((MI11*(U(i)-N11))-(MI12*N21)-(MI13*N31));
x5(i+1)=x5(i)+dt*((MI21*(U(i)-N11))-(MI22*N21)-(MI23*N31));
x6(i+1)=x6(i)+dt*((MI31*(U(i)-N11))-(MI32*N21)-(MI33*N31));
end
%fitness calculation(Using Energy of the system at final state of the experiment as a performance index. We want to minimize this) 
En1=(((1/2)*(m*((x4(701)^2)-(2*l1*(cos(x2(701)))*x4(701)*x5(701))+(11^2*x5(701)^2))))-m*g*l1*((cos(x2(701)))-1));
En2=(((1/2)*(m*((x4(701)^2)-(2*l2*(cos(x3(701)))*x4(701)*x6(701))+(l2^2*x6(701)^2))))-m*g*l2*((cos(x3(701)))-1));
En3=(1/2)*M*(x4(701)^2)+(1/2)*10*(x1(701)^2);
y=(((En1^2)+(En2^2)+(En3^2)));
