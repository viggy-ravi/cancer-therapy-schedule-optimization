%% Tirapazamine.m
%  Created by: Viggy Ravi
%  Last edited: 5/10/20

%% Tirapazamine
% COMMON PARAMETERS
k = 0.005;          % Proliferation rate (135 day doobling rate)
carcap = 8e8;      % Carrying Cap
days = 180;
dt = 1;

% BEVA SPECIFIC PARAMERTERS
beta = 0.1;        % Beva binding rate
VEGF = 2.59e-9;    % VEGF/cells      % Avg VEGF prod (Genentech)

% BEVA FACTS
ke_A = log(2)/20;  % beva elimination rate (half-life = 20 days, Genentech)

% VASCULATURE FACTS
cyto = 1200;      % mL/mg            % Cytotoxic effect (Imbs)
tau = 2;          % days             % delay parameter for Q (Imbs) 

% TIRA FACTS
CL = 898.848;       % L/day           % clearance (Senan)
Vdss = 39;          % L               % volume of distribution (Senan)
ke_C = CL/Vdss;     % 1/days          % rate of decay of Tira 
aC = 0.2;           % 1/day*t         % Tira internalization rate (Evans)

% TIRA INDUCED DEATH FACTS
O = 5;
dC = 4e9;             % carrying capacity reduction
S = 1;              % drug synergy

N  = zeros(1,days/dt);

A  = zeros(1,days/dt, 'double');
Af = zeros(1,days/dt, 'double');
Q  = zeros(1,days/dt);

h  = zeros(1,days/dt);
C  = zeros(1,days/dt);
Cf = zeros(1,days/dt);
pO2= 1.9*ones(1,days/dt);
oxy= zeros(1,days/dt);

t_A = [17 31 45]; 
t_C = t_A;      % time Tira is administered
beva = 10/1000;  % mg/g            % dosage of Beva 10 mg/g
tira = 0;       % mg/m^2          % dosage of Tira 330 mg/m^2 (Senan)

N(1) = 0.2*carcap;
Af(t_A/dt+1) = beva;
Cf(t_C/dt+1) = tira*1.7/1000;

idx=1;
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
   
   % Vasculature Quality - NOT SURE
   Q(t+tau/dt) = cyto*A(t-1);
   
   % pO2
   pO2(t) = 1.9 * (1+Q(t-1));
    
   % Tira concentration
   if t*dt - t_C(idx) > 0
      if Cf(t-1) - C(t-1) > 0
         Cf(t) = Cf(t-1) - C(t-1);
         C(t) = C(t-1) + dt*aC*Cf(t);
      else
         C(t) = C(t-1)*exp(-ke_C*(t*dt - t_C(idx)));   
      end

      % Tira induced death - NOT SURE
      oxy(t) = pO2(t)/O;
      if oxy(t) > 1
         oxy(t) = 1;
      end 
      h(t) = dC * C(t) * (1 - oxy(t));
   else
      C(t) = 0; 
      oxy(t) = pO2(t)/O;
      h(t) = 0;
   end

   % Tumor Growth   
   N(t) = N(t-1) + dt*(k*(1 - N(t-1)/(carcap-h(t)))*N(t-1) );
end

A = 1000.*A;
fprintf('Final num cells: %0.4f\n', N(days/dt)/reference)

%% Plots
uh = 1:1/dt:days/dt;
n = 3;

subplot(n,1,1)
plot(uh, N(uh))
xlabel('Days')
ylabel('Tumor Cells')
legend('Tumor Cells','Location', 'southeast')

% subplot(n,1,2)
% plot(uh, A(uh))
% xlabel('Days')
% ylabel('Beva Concentration ug/g')
% ylim([0 5])
% legend('Beva','Location', 'southeast')

subplot(n,1,2)
plot(uh, C(uh))
xlabel('Days')
ylabel('Tira Concentration mg/g')
ylim([0 0.5])
legend('Tira','Location', 'southeast')

subplot(n,1,3)
plot(uh, pO2(uh))
xlabel('Days')
ylabel('pO2 %')
legend('pO2','Location', 'southeast')

% subplot(n,1,4)
% plot(uh, h(uh))
% xlabel('Days')
% ylabel('h')
% legend('h','Location', 'southeast')
