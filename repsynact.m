load_syn=0;
load(sprintf('synapses%i.mat',load_syn),'synapses');

t = 0:0.00001:2.0;

figure(1)
%subplot(1,3,1)
for i=1:1:250
    
    i
    hold on
    plot(t,f(t,synapses(i)).*1e6);
    drawnow
end
xlabel('time [s]')
ylabel('influx [{\mu}mol/s.m2]')
%ylim([0 3.0e-8]*1e6)
title('Calcium Influx Pulses (all 250 synapses)');
set(gca,'FontSize',20)

% subplot(1,3,2)
% M = dlmread('aveCa.txt', ',', 1, 0);
% y = M(:,2);
% t = linspace(0,1,length(y));
% plot(t,y*1e6)
% xlabel('time [s]')
% ylabel('influx [{\mu}M]')
% ylim([0 9e-6]*1e6)
% title('Average Calcium for entire cell')
% set(gca,'FontSize',20)
% 
% subplot(1,3,3)
% y = M(:,4);
% plot(t,y)
% xlabel('time [s]')
% ylabel('influx [mV]')
% title('Average PM voltage for entire cell')
% set(gca,'FontSize',20)

set(gcf,'defaultAxesFontSize',20)

function y=f(t,s)
t0 = s.start_time;
tf = s.end_time;
m = -s.amp/(tf - t0);

y = t.*0;

%y((t>=t0)&(t<=tf)) = 1;
npulse = 2;
endtime = 1.0;

fr = endtime/npulse;

% twindow=[0.02 0.045];
% wend = twindow(2)+0.01;
% fr =wend + 1.2;

%fr = endtime/npulse;

for i=1:length(t)
    ti = mod(t(i),fr);
    
    if ((ti>=t0) & (ti<=tf))
        y(i) = m*(ti-t0) + s.amp;
    end    
end

end