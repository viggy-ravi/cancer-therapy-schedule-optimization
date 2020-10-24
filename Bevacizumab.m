%% Bevacizumab.m
%  Created by: Viggy Ravi
%  Last edited: 5/10/20

%% Bevacizumab
% COMMON PARAMETERS
k = 0.005;        % Proliferation rate (135 day doobling rate)
carcap = 8e8;    % Carrying Cap
days = 180;
dt = 1;

% BEVA SPECIFIC PARAMERTERS
beta = 0.1;        % Beva binding rate
VEGF = 2.59e-9;    % VEGF/cells      % Avg VEGF prod (Genentech)

% BEVA FACTS
dosage = 10/1000;  % mg/g            % dosage of Beva 10 mg/g
ke_A = log(2)/20;  % beva elimination rate (half-life = 20 days, Genentech)
t_A = [17 31 45]; 

% VASCULATURE FACTS
cyto = 1200;      % mL/mg            % Cytotoxic effect (Imbs)
tau = 2;          % days             % delay parameter for Q (Imbs) 

% TIRA (CHEMO) FACTS
gamma = 0.001;    % (g/mg)/day       % baseline effect of chemo (Imbs)
t_C = t_A + 4;

N  = zeros(1,days/dt);
A  = zeros(1,days/dt, 'double');
Af = zeros(1,days/dt, 'double');
Q  = zeros(1,days/dt);
pO2= ones(1,days/dt);

N(1) = 0.2*carcap;
Af(t_A/dt+1) = dosage;

idx = 1;
for t = 2:(days/dt)
   
   % find dosage
   for i = 1:length(t_A)
      if t==t_A(i)/dt 
         idx=i;
      end
   end
   % Beva concentration
   if Af(t-1) - A(t-1) > 0
      Af(t) = Af(t-1) - A(t-1);
      A(t) = A(t-1) + dt*(beta*Af(t)*(N(t-1)*VEGF - A(t-1)));
   else
      A(t) = A(t-1)*exp(-ke_A*((t*dt) - t_A(idx)));   
   end
   
   % Vasculature Quality
   Q(t+tau/dt) = cyto*A(t-1);
   
   pO2(t) = 1.9 * (1+Q(t-1));
   % Tumor Growth   
   N(t) = N(t-1) + dt*(k*(1 - N(t-1)/carcap)*N(t-1) - gamma*Q(t-1)*N(t-1));
end

A = 1000.*A;

%% Plot
uh = 60;

subplot(4,1,1)
p = plot(1:uh, N(1:uh));
p(1).LineWidth = 3;
xlabel('Days')
ylabel('Tumor Cells')
% legend('Tumor Cells','Location', 'southeast')
title('Multiple Dosages of Beva - Days 17,31,45')

subplot(4,1,2)
p = plot(1:uh, A(1:uh), 'b');
p(1).LineWidth = 3;
xlabel('Days')
xticks(0:5:uh)
ylabel('Beva Conc. (ug/g)')
% legend('Beva','Location', 'southwest')
title('Bevacizumab Concentration')

subplot(4,1,3)
p = plot(1:uh, Q(1:uh),'r');
p(1).LineWidth = 3;
xlabel('Days')
xticks(0:5:uh)
ylabel('Vasc. Qual.')
% legend('Vasc. Qual.','Location', 'southwest')
title('Vasculature Quality')

subplot(4,1,4)
p = plot(1:uh, pO2(1:uh),'g');
p(1).LineWidth = 3;
xlabel('Days')
ylabel('pO2 (%)')
% legend('pO2 %','Location', 'southwest')
title('Oxygenation')