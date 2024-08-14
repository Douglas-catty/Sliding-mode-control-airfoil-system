clear
clc

h=1e-3;
d=0.07;
Q=8;
k1=20;

A21=(-1/1.75)*(2-0.1*Q-4*d*Q);
A22=(-1/1.75)*(0.4);
A23=(-1/1.75)*(-0.2);
A24=(-1/1.75)*(-0.1);
A25=(-80/1.75);

A41=(-1/1.75)*(-0.5+0.2*Q+d*Q);
A42=(-1/1.75)*(-0.1);
A43=(-1/1.75)*(0.4);
A44=(-1/1.75)*(0.2);
A45=20/1.75;

%x10=0.5;x20=0.5;x30=0.5;x40=0.5;
x10=0.1331;x20=0.1265;x30=0.4508;x40=0.1198;
n=2e5;
x=zeros(n,4);
x(1,:)=[x10,x20,x30,x40];

% Airfoil system
f1=@(x,y,z,m) y;
f2=@(x,y,z,m) A21*x+A22*y+A23*z+A24*m+A25*x^3;
f3=@(x,y,z,m) m;
f4=@(x,y,z,m) A41*x+A42*y+A43*z+A44*m+A45*x^3;

% Noise parameters
 % sigmaB1=0.02;
sigmaB2=0.1;
 % sigmaB3=0.02;
sigmaB4=0.1;

alpha1Levy=1.75;
alpha2Levy=1.75;
alpha3Levy=1.75;
alpha4Levy=1.75;
beta1Levy=0;
beta2Levy=0;
beta3Levy=0;
beta4Levy=0;

sigma1Levy=0.01;  % Levy noise intensity
sigma2Levy=0.01;
sigma3Levy=0.01;
sigma4Levy=0.01;


% Discrete control

tend=h*n;
tline=linspace(h,tend,tend/h);

omiga=1;
r=0.3;

%%%%%%%%%%%%%%%%%%%%%%%%%
Example=4;

if Example==1

Target1=@(t) r*cos(omiga*t);
Target2=@(t) r*sin(omiga*t);
Target1D1=@(t) -r*omiga*sin(omiga*t);
Target2D1=@(t)  r*omiga*cos(omiga*t);
Target1D2=@(t) -r*omiga*omiga*cos(omiga*t);
Target2D2=@(t) -r*omiga*omiga*sin(omiga*t);

elseif Example==2

Target1=@(t) 0.05*(cos(t)*(exp(sin(8*t)) + 1));
Target2=@(t) 0.05*(sin(t)*(exp(sin(8*t)) + 1));
Target1D1=@(t) 0.05*(8*cos(8*t)*exp(sin(8*t))*cos(t) - sin(t)*(exp(sin(8*t)) + 1));
Target2D1=@(t) 0.05*(cos(t)*(exp(sin(8*t)) + 1) + 8*cos(8*t)*exp(sin(8*t))*sin(t)); 
Target1D2=@(t) 0.05*(64*cos(8*t)^2*exp(sin(8*t))*cos(t) - 16*cos(8*t)*exp(sin(8*t))*sin(t) - 64*sin(8*t)*exp(sin(8*t))*cos(t) - cos(t)*(exp(sin(8*t)) + 1));
Target2D2=@(t) 0.05*(16*cos(8*t)*exp(sin(8*t))*cos(t) - sin(t)*(exp(sin(8*t)) + 1) - 64*sin(8*t)*exp(sin(8*t))*sin(t) + 64*cos(8*t)^2*exp(sin(8*t))*sin(t));

elseif Example==3

Target1=@(t) 0.01*(16*(sin(t)^3));
Target2=@(t) 0.01*(13*cos(t)-5*cos(2*t)-2*cos(3*t)-cos(4*t));
Target1D1=@(t) 0.01*(48*cos(t)*sin(t)^2);
Target2D1=@(t) 0.01*(10*sin(2*t) + 6*sin(3*t) + 4*sin(4*t) - 13*sin(t)); 
Target1D2=@(t) 0.01*(96*cos(t)^2*sin(t) - 48*sin(t)^3);
Target2D2=@(t) 0.01*(20*cos(2*t) + 18*cos(3*t) + 16*cos(4*t) - 13*cos(t));

elseif Example==4

Target1=@(t) 0.1*(-sin(t)-sin(t/2));
Target2=@(t) 0.1*(-cos(t)+cos(t/2));
Target1D1=@(t) 0.1*(- cos(t/2)/2 - cos(t));
Target2D1=@(t) 0.1*(sin(t) - sin(t/2)/2); 
Target1D2=@(t) 0.1*(sin(t/2)/4 + sin(t));
Target2D2=@(t) 0.1*(cos(t) - cos(t/2)/4);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%


k1=10;  % Designed exponential convergence rate
k2=10;
k3=10;
k4=10;

for k=2:n

    [k,n]

xnow=x(k-1,1);
ynow=x(k-1,2);
znow=x(k-1,3);
mnow=x(k-1,4);
tnow=tline(1,k-1);

%%% RK4
K1=f1(xnow,ynow,znow,mnow);
K2=f1(xnow+h*K1/2,ynow+h*K1/2,znow+h*K1/2,mnow+h*K1/2);
K3=f1(xnow+h*K2/2,ynow+h*K2/2,znow+h*K2/2,mnow+h*K2/2);
K4=f1(xnow+h*K3,ynow+h*K3,znow+h*K3,mnow+h*K3);
L1=f2(xnow,ynow,znow,mnow);
L2=f2(xnow+h*L1/2,ynow+h*L1/2,znow+h*L1/2,mnow+h*L1/2);
L3=f2(xnow+h*L2/2,ynow+h*L2/2,znow+h*L2/2,mnow+h*L2/2);
L4=f2(xnow+h*L3,ynow+h*L3,znow+h*L3,mnow+h*L3);
M1=f3(xnow,ynow,znow,mnow);
M2=f3(xnow+h*M1/2,ynow+h*M1/2,znow+h*M1/2,mnow+h*M1/2);
M3=f3(xnow+h*M2/2,ynow+h*M2/2,znow+h*M2/2,mnow+h*M2/2);
M4=f3(xnow+h*M3,ynow+h*M3,znow+h*M3,mnow+h*M3);
N1=f4(xnow,ynow,znow,mnow);
N2=f4(xnow+h*N1/2,ynow+h*N1/2,znow+h*N1/2,mnow+h*N1/2);
N3=f4(xnow+h*N2/2,ynow+h*N2/2,znow+h*N2/2,mnow+h*N2/2);
N4=f4(xnow+h*N3,ynow+h*N3,znow+h*N3,mnow+h*N3);

dx=(K1+2*K2+2*K3+K4)*h/6;
dy=(L1+2*L2+2*L3+L4)*h/6;
dz=(M1+2*M2+2*M3+M4)*h/6;
dm=(N1+2*N2+2*N3+N4)*h/6;

% Noise
Bh=sqrt(h)*randn(4,1);

Mx=stblrnd(alpha1Levy,beta1Levy,h^(1/alpha1Levy),0,1,1);
My=stblrnd(alpha2Levy,beta2Levy,h^(1/alpha2Levy),0,1,1);
Mz=stblrnd(alpha3Levy,beta3Levy,h^(1/alpha3Levy),0,1,1);
Mm=stblrnd(alpha4Levy,beta4Levy,h^(1/alpha4Levy),0,1,1);
Levyh=[sigma1Levy*Mx;sigma2Levy*My;sigma3Levy*Mz;sigma4Levy*Mm];

    %x(k,1)=x(k-1,1)+dx+sigmaB1*Bh(1,1);
    x(k,1)=x(k-1,1)+dx;
    x(k,2)=x(k-1,2)+dy+sigmaB2*Bh(2,1)+Levyh(2,1);
    %x(k,2)=x(k-1,2)+dy-h*epsilon*sign(s1(x(k-1,1),x(k-1,2),x(k-1,3),x(k-1,4)))+sigmaB2*Bh(2,1)+Levyh(2,1);
    
    
    x(k,3)=x(k-1,3)+dz;
    x(k,4)=x(k-1,4)+dm+sigmaB4*Bh(4,1)+Levyh(4,1);
    
    u2=k2*(k1*(Target1(tnow)-xnow)+Target1D1(tnow)-ynow)+Target1(tnow)-xnow+k1*(Target1D1(tnow)-ynow)+Target1D2(tnow)...
        -(A21*xnow+A22*ynow+A23*znow+A24*mnow+A25*xnow^3);

    u4=k4*(k3*(Target2(tnow)-znow)+Target2D1(tnow)-mnow)+Target2(tnow)-znow+k3*(Target2D1(tnow)-mnow)+Target2D2(tnow)...
        -(A41*xnow+A42*ynow+A43*znow+A44*mnow+A45*xnow^3);


    x(k,2)=x(k,2)+u2*h;
    x(k,4)=x(k,4)+u4*h;
    
end


subplot(3,1,1)
set(gcf,'color','white');
plot(x(:,1),x(:,3),'b','Linewidth',0.5);
hold on

load PeriodicTrajectory.mat

plot(PeriodicTrajectory(:,1),PeriodicTrajectory(:,3),'r','Linewidth',1.5);
hold on

xlabel('x_1');
ylabel('x_3');
legend('Controlled trajectory','Flutter periodic orbit');

subplot(3,1,2)

if Example==1

plot(tline,r*cos(omiga*tline),'b','Linewidth',0.5);
hold on
plot(tline,x(:,1),'r--','Linewidth',0.5);
hold on

elseif Example==2

plot(tline,0.05*(cos(tline).*(exp(sin(8*tline)) + 1)),'b','Linewidth',0.5);
hold on
plot(tline,x(:,1),'r--','Linewidth',0.5);
hold on

elseif Example==3

plot(tline,0.01*(16*(sin(tline).^3)),'b','Linewidth',0.5);
hold on
plot(tline,x(:,1),'r--','Linewidth',0.5);
hold on

elseif Example==4

plot(tline,0.1*(-sin(tline)-sin(tline/2)),'b','Linewidth',0.5);
hold on
plot(tline,x(:,1),'r--','Linewidth',0.5);
hold on

end

xlabel('t');
ylabel('x_1');
legend('Target','Controlled trajectory');


subplot(3,1,3)

if Example==1

plot(tline,r*sin(omiga*tline),'b','Linewidth',0.5);
hold on
plot(tline,x(:,3),'r--','Linewidth',0.5);
hold on

elseif Example==2

plot(tline,0.05*(sin(tline).*(exp(sin(8*tline)) + 1)),'b','Linewidth',0.5);
hold on
plot(tline,x(:,3),'r--','Linewidth',0.5);
hold on

elseif Example==3

plot(tline,0.01*(13*cos(tline)-5*cos(2*tline)-2*cos(3*tline)-cos(4*tline)),'b','Linewidth',0.5);
hold on
plot(tline,x(:,3),'r--','Linewidth',0.5);
hold on

elseif Example==4

plot(tline,0.1*(-cos(tline)+cos(tline/2)),'b','Linewidth',0.5);
hold on
plot(tline,x(:,3),'r--','Linewidth',0.5);
hold on

end

xlabel('t');
ylabel('x_3');
legend('Target','Controlled trajectory');



