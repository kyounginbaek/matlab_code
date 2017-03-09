% +++++++++++++++++++
% 
% Matlab codes for HH model.
% J Rinzel 2006  (in part, from D Tranchina)
% +++++++++++++++++++

%F=farads
%V=mV
%A=amps
%F=Farads
C=1.0;       %micro F/cm^2
gbar_Na=120; %(micro A/mV)/cm^2
gbar_K=36;   %(micro A/mV)/cm^2
gbar_L=0.3;   %(micro A/mV)/cm^2
E_Na=45;     %mV
E_K=-82;     %mV
E_L=-59;     %mv
%V=Y(1);
%m=Y(2);
%h=Y(3);
%n=Y(4);
Y0=[-69.89825095976047, 0.05357073962900, 0.59257698454360, 0.31924676125340];
%Y0(1)=Y0(1)+10;
t=0:0.01:60;
options=odeset('RelTol',1.e-8);
[T, Y]=ode45(@dydt_HH,t,Y0,options, C,gbar_Na,gbar_K,gbar_L,E_Na,E_K,E_L,...
             'Ioft');
%JR figure``(1)
figure(1);
plot(T,Y(:,1)) % time course of V(t)
xlabel('Time (ms)')
ylabel('V (mV)')

figure(2);
plot(T,(Y(:,1)+70)./100,T,Y(:,2),T,Y(:,3),T,Y(:,4)); % time courses of V (scaled),m,h,and n
legend('scaled V(t)','m(t)','h(t)','n(t)');
xlabel('Time'); ylabel('V and gating vars');

figure(3); clf;
hold on;
plot(Y(:,1),Y(:,2)); % V-m phase plane projection
Vpts=(-100.001:.5:40);
minfpts=m_and_tau_m(Vpts);
plot(Vpts,minfpts,'black');
xlabel('V (mV)'); ylabel('m');

figure(4); clf;
hold on;
plot(Y(:,1),Y(:,4)); % V-n phase plane projection
ninfpts=n_and_tau_n(Vpts);
plot(Vpts,ninfpts,'black');
xlabel('V (mV)'); ylabel('n');

figure(5);
plot(Y(:,4),Y(:,3)); % n-h phase plane projection
xlabel('n'); ylabel('h');

