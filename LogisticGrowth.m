% Logistic Growth - NSCLC

k = 0.005;        % Proliferation rate (doubling rate = 135 days)
carcap = 8e8;     % Carrying Cap

dt = 1;
days = 180;       %0.5 years


N = zeros(1,days/dt);
N(1) = 0.20*carcap;

for t = 2:(days/dt)   
   % Tumor Growth
   N(t) = N(t-1) + dt*(k*(1 - N(t-1)/carcap)*N(t-1));
end

%% Plots

uh = 1:1/dt:days/dt;
n = 1;

subplot(n,1,1)
plot(uh, N(uh))
xlabel('Days')
ylabel('Tumor Cells')
ylim([0,carcap])
legend('Tumor Cells','Location', 'southeast')