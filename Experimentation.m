%% Experimentation
clc; clear all;
load('beva_tira')
trials = 70;

N = zeros(trials,days/dt);
A = zeros(trials,days/dt);
C = zeros(trials,days/dt);
r = zeros(trials);

% Logistic Growth
n=1;
[N(n,:),A(n,:),C(n,:)] = sequence_optimization(0,0,[0],[0]);
reference = N(1,days/dt);
r(n) = N(n,days/dt)/reference;
fprintf('Control: \t\t%0.4f\n', N(n,days/dt)/reference)

% Individual (Therapies Alone)
n=2;
dosage_times = [0];
delays = [0];

[N(n,:),A(n,:),C(n,:)] = sequence_optimization(10,0,dosage_times,delays);
r(n) = N(n,days/dt)/reference;
fprintf('Beva Only: \t\t%0.4f\n', N(n,days/dt)/reference)

n=3;
[N(n,:)] = sequence_optimization(0,330,dosage_times,delays);
r(n) = N(n,days/dt)/reference;
fprintf('Tira Only: \t\t%0.4f\n', N(n,days/dt)/reference)

% Beva-Tira 1 dose
n=4:12;
dosage_times = [17];
delays = -4:4;
[N(n,:),A(n,:),C(n,:)] = sequence_optimization(10,330,dosage_times,delays);
r(n) = N(n,days/dt)/reference;
for i=n
   fprintf('%iDose (%i): \t%0.4f\n', length(dosage_times), delays(i-(n(1)-1)), N(i,days/dt)/reference) 
end

% Beva-Tira 2 doses
n=13:21;
dosage_times = [17 31];
delays = -4:4;
[N(n,:),A(n,:),C(n,:)] = sequence_optimization(10,330,dosage_times,delays);
r(n) = N(n,days/dt)/reference;
for i=n
   fprintf('%iDose (%i): \t%0.4f\n', length(dosage_times), delays(i-(n(1)-1)), N(i,days/dt)/reference) 
end

% Beva-Tira 3 doses
n=22:30;
dosage_times = [17 31 45];
delays = -4:4;
[N(n,:),A(n,:),C(n,:)] = sequence_optimization(10,330,dosage_times,delays);
r(n) = N(n,days/dt)/reference;
for i=n
   fprintf('%iDose (%i): \t%0.4f\n', length(dosage_times), delays(i-(n(1)-1)), N(i,days/dt)/reference) 
end

% Beva-Tira 4 doses
n=31:39;
dosage_times = [17 31 45 59];
delays = -4:4;
[N(n,:),A(n,:),C(n,:)] = sequence_optimization(10,330,dosage_times,delays);
r(n) = N(n,days/dt)/reference;
for i=n
   fprintf('%iDose (%i): \t%0.4f\n', length(dosage_times), delays(i-(n(1)-1)), N(i,days/dt)/reference) 
end

% Beva-Tira 5 doses
n=40:48;
dosage_times = [17 31 45 59 73];
delays = -4:4;
[N(n,:),A(n,:),C(n,:)] = sequence_optimization(10,330,dosage_times,delays);
r(n) = N(n,days/dt)/reference;
for i=n
   fprintf('%iDose (%i): \t%0.4f\n', length(dosage_times), delays(i-(n(1)-1)), N(i,days/dt)/reference) 
end

% Beva-Tira 6 doses
n=49:57;
dosage_times = [17 31 45 59 73 87];
delays = -4:4;
[N(n,:),A(n,:),C(n,:)] = sequence_optimization(10,330,dosage_times,delays);
r(n) = N(n,days/dt)/reference;
for i=n
   fprintf('%iDose (%i): \t%0.4f\n', length(dosage_times), delays(i-(n(1)-1)), N(i,days/dt)/reference) 
end

%% Plots
t = 1:dt:60;

figure(1)
p = plot(t,N(1,t),'k', t,N(2,t), t,N(3,t));
p(1).LineWidth = 2;
title('Individually - Administered on Day 0')
xlabel('Days')
ylabel('Num. Tumor Cells')
% ylim([0,2e8])
legend('Control','Beva Only','Tira Only','Location','southeast')

figure(9)
p = plot(t,N(1,t),'k');
p(1).LineWidth = 2;
title('Control')
xlabel('Days')
ylabel('Num. Tumor Cells')
% ylim([0,2e8])
legend('Control','Location','southeast')

figure(2)
plot(t,N(1,t), t,N(4:12,t))
title('1 Dose')
xlabel('Days')
ylabel('Num. Tumor Cells')
legend('Control','-4','-3','-2','-1','0','1','2','3','4','Location','southeast')

figure(3)
plot(t,N(1,t), t,N(13:21,t))
title('2 Doses')
xlabel('Days')
ylabel('Num. Tumor Cells')
legend('Control','-4','-3','-2','-1','0','1','2','3','4','Location','southeast')

figure(4)
plot(t,N(1,t), t,N(22:30,t))
title('3 Doses')
xlabel('Days')
ylabel('Num. Tumor Cells')
legend('Control','-4','-3','-2','-1','0','1','2','3','4','Location','southeast')

figure(5)
plot(t,N(1,t), t,N(31:39,t))
title('4 Doses')
xlabel('Days')
ylabel('Num. Tumor Cells')
legend('Control','-4','-3','-2','-1','0','1','2','3','4','Location','southeast')

figure(6)
plot(t,N(1,t), t,N(40:48,t))
title('5 Doses')
xlabel('Days')
ylabel('Num. Tumor Cells')
legend('Control','-4','-3','-2','-1','0','1','2','3','4','Location','southeast')

figure(7)
p = plot(t,N(1,t),'k', t,N(49,t),'--', t,N(52,t), t,N(54,t), t,N(57,t),'--');
p(1).LineWidth = 2;
p(4).LineWidth = 2;
title('6 Doses')
xlabel('Days')
xticks(0:5:days/dt)
ylabel('Num. Tumor Cells')
legend('Control','Tira-Beva 4Days','Tira-Beva 1Day',...
       'Beva-Tira 1Day','Beva-Tira 4Days','Location','northwest')

figure(8)
best = 9:9:54;
plot(t,N(1,t), t,N(best,t))
title('Beva-Tira:Delay +1')
xlabel('Days')
ylabel('Num. Tumor Cells')
legend('Control','1Dose','2Dose','3Dose','4Dose','5Dose','6Dose','Location','southeast')
