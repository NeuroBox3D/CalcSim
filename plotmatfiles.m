load('cc.mat');
load('ce.mat');
load('bb.mat');
load('bbe.mat');
load('vv.mat');
load('o1.mat');
load('o2.mat');

nT = length(mat1(1,:));
t = linspace(0,1,nT);

figure(1)
xn = 250;

subplot(1,7,1);
plot(t,mat1(xn,:)*1e6)
title('Cyt [Ca2+]')
xlabel('time [s]')
ylabel(' {\mu}M')

subplot(1,7,2);
plot(t,mat2(xn,:)*1e6)
title('Cyt [CalB]')
xlabel('time [s]')
ylabel(' {\mu}M')

subplot(1,7,3);
plot(t,mat3(xn,:)*1e6)
title('ER [Ca2+]')
xlabel('time [s]')
ylabel(' {\mu}M')

subplot(1,7,4);
plot(t,mat4(xn,:)*1e6)
title('ER Buff')
xlabel('time [s]')
ylabel(' {\mu}M')

subplot(1,7,5);
plot(t,mat5(xn,:))
title('PM Voltage')
xlabel('time [s]')
ylabel('mV')

subplot(1,7,6);
plot(t,mat6(xn,:))
title('O1 State')
ylim([0 1])
xlabel('time [s]')
ylabel('[]')

subplot(1,7,7);
plot(t,mat7(xn,:))
title('O2 State')
ylim([0 1])
xlabel('time [s]')
ylabel('[]')