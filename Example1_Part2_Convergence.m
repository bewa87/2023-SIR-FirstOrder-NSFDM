%
% Example 1: SIR Model - Convergence of NSFDM
%

% Step 0: Definition Of Problem Parameters

T       = 70;                                     % Final Simulation T
h       = 0.0001;                                 % Time Step Size For RK4
t       = 0:h:T;                                  % Time Vector For RK4
h_nsfdm = [1; 0.5; 0.25; 0.125; 0.0625; 0.03125]; % Time Step Sizes For NSFDM
t1      = 0:h_nsfdm(1):T;                         % Time Vector For 1. NSFDM
t2      = 0:h_nsfdm(2):T;                         % Time Vector For 2. NSFDM
t3      = 0:h_nsfdm(3):T;                         % Time Vector For 3. NSFDM
t4      = 0:h_nsfdm(4):T;                         % Time Vector For 4. NSFDM
t5      = 0:h_nsfdm(5):T;                         % Time Vector For 5. NSFDM
t6      = 0:h_nsfdm(6):T;                         % Time Vector For 6. NSFDM
alpha   = 0.5;                                    % Time-Varying Transmission Rate alpha
beta    = 0.1;                                    % Time-Varying Recovery Rate beta
S_1     = 5000;                                   % Initial Condition For Susceptible People
I_1     = 5000;                                   % Initial Condition For Infected People
R_1     = 0;                                      % Initial Condition For Recovered People
N       = S_1 + I_1 + R_1;                        % Constant Total Population Size N

% Step 1: Initialization Of Solution Vectors

% Step 1.1: Initialization For RK4

S_RK4    = zeros(length(t),1);
I_RK4    = zeros(length(t),1);
R_RK4    = zeros(length(t),1);
S_RK4(1) = S_1;
I_RK4(1) = I_1;
R_RK4(1) = R_1;

% Step 1.2: Initialization For 1. NSFDM

S_NSFDM1    = zeros(length(t1),1);
I_NSFDM1    = zeros(length(t1),1);
R_NSFDM1    = zeros(length(t1),1);
S_NSFDM1(1) = S_1;
I_NSFDM1(1) = I_1;
R_NSFDM1(1) = R_1;

% Step 1.3: Initialization For 2. NSFDM

S_NSFDM2    = zeros(length(t2),1);
I_NSFDM2    = zeros(length(t2),1);
R_NSFDM2    = zeros(length(t2),1);
S_NSFDM2(1) = S_1;
I_NSFDM2(1) = I_1;
R_NSFDM2(1) = R_1;

% Step 1.4: Initialization For 3. NSFDM

S_NSFDM3    = zeros(length(t3),1);
I_NSFDM3    = zeros(length(t3),1);
R_NSFDM3    = zeros(length(t3),1);
S_NSFDM3(1) = S_1;
I_NSFDM3(1) = I_1;
R_NSFDM3(1) = R_1;

% Step 1.5: Initialization For 4. NSFDM

S_NSFDM4    = zeros(length(t4),1);
I_NSFDM4    = zeros(length(t4),1);
R_NSFDM4    = zeros(length(t4),1);
S_NSFDM4(1) = S_1;
I_NSFDM4(1) = I_1;
R_NSFDM4(1) = R_1;

% Step 1.6: Initialization For 5. NSFDM

S_NSFDM5    = zeros(length(t5),1);
I_NSFDM5    = zeros(length(t5),1);
R_NSFDM5    = zeros(length(t5),1);
S_NSFDM5(1) = S_1;
I_NSFDM5(1) = I_1;
R_NSFDM5(1) = R_1;

% Step 1.7: Initialization For 3. NSFDM

S_NSFDM6    = zeros(length(t6),1);
I_NSFDM6    = zeros(length(t6),1);
R_NSFDM6    = zeros(length(t6),1);
S_NSFDM6(1) = S_1;
I_NSFDM6(1) = I_1;
R_NSFDM6(1) = R_1;

% Step 2: Loops Over All Time Points For All Different h

% Step 2.1: Loop For RK4

for j = 1:1:length(t)-1
  k1 = [-alpha*I_RK4(j)*S_RK4(j)/N; alpha*I_RK4(j)*S_RK4(j)/N-beta*I_RK4(j); beta*I_RK4(j)];
  k2 = [-alpha*((I_RK4(j)+h/2*k1(2))*(S_RK4(j)+h/2*k1(1)))/N; ...
         alpha*((I_RK4(j)+h/2*k1(2))*(S_RK4(j)+h/2*k1(1)))/N-beta*(I_RK4(j)+h/2*k1(2)); ...
         beta*(I_RK4(j)+h/2*k1(2))];
  k3 = [-alpha*((I_RK4(j)+h/2*k2(2))*(S_RK4(j)+h/2*k2(1)))/N; ...
         alpha*((I_RK4(j)+h/2*k2(2))*(S_RK4(j)+h/2*k2(1)))/N-beta*(I_RK4(j)+h/2*k2(2)); ...
         beta*(I_RK4(j)+h/2*k2(2))];
  k4 = [-alpha*((I_RK4(j)+h*k3(2))*(S_RK4(j)+h*k3(1)))/N; ...
         alpha*((I_RK4(j)+h*k3(2))*(S_RK4(j)+h*k3(1)))/N-beta*(I_RK4(j)+h*k3(2)); ...
         beta*(I_RK4(j)+h*k3(2))];
  S_RK4(j+1) = S_RK4(j) + h/6*(k1(1)+2*k2(1)+2*k3(1)+k4(1));
  I_RK4(j+1) = I_RK4(j) + h/6*(k2(2)+2*k2(2)+2*k3(2)+k4(2));
  R_RK4(j+1) = R_RK4(j) + h/6*(k2(3)+2*k3(3)+2*k3(3)+k4(3));
endfor

% Step 2.2: Loop For 1. NSFDM

for j = 1:1:length(t1)-1
  S_NSFDM1(j+1) = S_NSFDM1(j)/(1+h_nsfdm(1)*alpha*I_NSFDM1(j)/N);
  I_NSFDM1(j+1) = (I_NSFDM1(j)+h_nsfdm(1)*alpha*S_NSFDM1(j+1)*I_NSFDM1(j)/N)/(1+h_nsfdm(1)*beta);
  R_NSFDM1(j+1) = R_NSFDM1(j)+h_nsfdm(1)*beta*I_NSFDM1(j+1);
endfor

% Step 2.3: Loop For 2. NSFDM

for j = 1:1:length(t2)-1
  S_NSFDM2(j+1) = S_NSFDM2(j)/(1+h_nsfdm(2)*alpha*I_NSFDM2(j)/N);
  I_NSFDM2(j+1) = (I_NSFDM2(j)+h_nsfdm(2)*alpha*S_NSFDM2(j+1)*I_NSFDM2(j)/N)/(1+h_nsfdm(2)*beta);
  R_NSFDM2(j+1) = R_NSFDM2(j)+h_nsfdm(2)*beta*I_NSFDM2(j+1);
endfor

% Step 2.4: Loop For 3. NSFDM

for j = 1:1:length(t3)-1
  S_NSFDM3(j+1) = S_NSFDM3(j)/(1+h_nsfdm(3)*alpha*I_NSFDM3(j)/N);
  I_NSFDM3(j+1) = (I_NSFDM3(j)+h_nsfdm(3)*alpha*S_NSFDM3(j+1)*I_NSFDM3(j)/N)/(1+h_nsfdm(3)*beta);
  R_NSFDM3(j+1) = R_NSFDM3(j)+h_nsfdm(3)*beta*I_NSFDM3(j+1);
endfor

% Step 2.5: Loop For 4. NSFDM

for j = 1:1:length(t4)-1
  S_NSFDM4(j+1) = S_NSFDM4(j)/(1+h_nsfdm(4)*alpha*I_NSFDM4(j)/N);
  I_NSFDM4(j+1) = (I_NSFDM4(j)+h_nsfdm(4)*alpha*S_NSFDM4(j+1)*I_NSFDM4(j)/N)/(1+h_nsfdm(4)*beta);
  R_NSFDM4(j+1) = R_NSFDM4(j)+h_nsfdm(4)*beta*I_NSFDM4(j+1);
endfor

% Step 2.6: Loop For 5. NSFDM

for j = 1:1:length(t5)-1
  S_NSFDM5(j+1) = S_NSFDM5(j)/(1+h_nsfdm(5)*alpha*I_NSFDM5(j)/N);
  I_NSFDM5(j+1) = (I_NSFDM5(j)+h_nsfdm(5)*alpha*S_NSFDM5(j+1)*I_NSFDM5(j)/N)/(1+h_nsfdm(5)*beta);
  R_NSFDM5(j+1) = R_NSFDM5(j)+h_nsfdm(5)*beta*I_NSFDM5(j+1);
endfor

% Step 2.7: Loop For 6. NSFDM

for j = 1:1:length(t6)-1
  S_NSFDM6(j+1) = S_NSFDM6(j)/(1+h_nsfdm(6)*alpha*I_NSFDM6(j)/N);
  I_NSFDM6(j+1) = (I_NSFDM6(j)+h_nsfdm(6)*alpha*S_NSFDM6(j+1)*I_NSFDM6(j)/N)/(1+h_nsfdm(6)*beta);
  R_NSFDM6(j+1) = R_NSFDM6(j)+h_nsfdm(6)*beta*I_NSFDM6(j+1);
endfor

% Step 3: Summary Of Values For Errors

h_sizes = [0.03125 0.0625 0.125 0.25 0.5 1];
errors_S = [abs(S_NSFDM6(end)-S_RK4(end)) ...
            abs(S_NSFDM5(end)-S_RK4(end)) ... 
            abs(S_NSFDM4(end)-S_RK4(end)) ... 
            abs(S_NSFDM3(end)-S_RK4(end)) ... 
            abs(S_NSFDM2(end)-S_RK4(end)) ...
            abs(S_NSFDM1(end)-S_RK4(end))];
errors_I = [abs(I_NSFDM6(end)-I_RK4(end)) ...
            abs(I_NSFDM5(end)-I_RK4(end)) ... 
            abs(I_NSFDM4(end)-I_RK4(end)) ... 
            abs(I_NSFDM3(end)-I_RK4(end)) ... 
            abs(I_NSFDM2(end)-I_RK4(end)) ...
            abs(I_NSFDM1(end)-I_RK4(end))];
errors_R = [abs(R_NSFDM6(end)-R_RK4(end)) ...
            abs(R_NSFDM5(end)-R_RK4(end)) ... 
            abs(R_NSFDM4(end)-R_RK4(end)) ... 
            abs(R_NSFDM3(end)-R_RK4(end)) ... 
            abs(R_NSFDM2(end)-R_RK4(end)) ...
            abs(R_NSFDM1(end)-R_RK4(end))];