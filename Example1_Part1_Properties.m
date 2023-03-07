%
% Example 1: SIR Model - Comparison of NSFDM, EE, RK2
%

% Step 0: Definition Of Problem Parameters

T     = 70;              % Final Simulation T
h1    = 7.5;             % Time Step Size h1
t1    = 0:h1:T;          % Time Vector t1
h2    = 5;               % Time Step Size h2
t2    = 0:h2:T;          % Time Vector t2
h3    = 0.1;             % Time Step Size h3
t3    = 0:h3:T;          % Time Vector t3
alpha = 0.5;             % Time-Varying Transmission Rate alpha
beta  = 0.1;             % Time-Varying Recovery Rate beta
S_1   = 5000;            % Initial Condition For Susceptible People
I_1   = 5000;            % Initial Condition For Infected People
R_1   = 0;               % Initial Condition For Recovered People
N     = S_1 + I_1 + R_1; % Constant Total Population Size N

% Step 1: Initialization Of Solution Vectors

% Step 1.1: Initialization For NSFDM

S_NSFDM1    = zeros(length(t1),1);
I_NSFDM1    = zeros(length(t1),1);
R_NSFDM1    = zeros(length(t1),1);
S_NSFDM1(1) = S_1;
I_NSFDM1(1) = I_1;
R_NSFDM1(1) = R_1;
S_NSFDM2    = zeros(length(t2),1);
I_NSFDM2    = zeros(length(t2),1);
R_NSFDM2    = zeros(length(t2),1);
S_NSFDM2(1) = S_1;
I_NSFDM2(1) = I_1;
R_NSFDM2(1) = R_1;
S_NSFDM3    = zeros(length(t3),1);
I_NSFDM3    = zeros(length(t3),1);
R_NSFDM3    = zeros(length(t3),1);
S_NSFDM3(1) = S_1;
I_NSFDM3(1) = I_1;
R_NSFDM3(1) = R_1;

% Step 1.2: Initialization For EE

S_EE1    = zeros(length(t1),1);
I_EE1    = zeros(length(t1),1);
R_EE1    = zeros(length(t1),1);
S_EE1(1) = S_1;
I_EE1(1) = I_1;
R_EE1(1) = R_1;
S_EE2    = zeros(length(t2),1);
I_EE2    = zeros(length(t2),1);
R_EE2    = zeros(length(t2),1);
S_EE2(1) = S_1;
I_EE2(1) = I_1;
R_EE2(1) = R_1;
S_EE3    = zeros(length(t3),1);
I_EE3    = zeros(length(t3),1);
R_EE3    = zeros(length(t3),1);
S_EE3(1) = S_1;
I_EE3(1) = I_1;
R_EE3(1) = R_1;

% Step 1.3: Initialization For RK2

S_RK2_1    = zeros(length(t1),1);
I_RK2_1    = zeros(length(t1),1);
R_RK2_1    = zeros(length(t1),1);
S_RK2_1(1) = S_1;
I_RK2_1(1) = I_1;
R_RK2_1(1) = R_1;
S_RK2_2    = zeros(length(t2),1);
I_RK2_2    = zeros(length(t2),1);
R_RK2_2    = zeros(length(t2),1);
S_RK2_2(1) = S_1;
I_RK2_2(1) = I_1;
R_RK2_2(1) = R_1;
S_RK2_3    = zeros(length(t3),1);
I_RK2_3    = zeros(length(t3),1);
R_RK2_3    = zeros(length(t3),1);
S_RK2_3(1) = S_1;
I_RK2_3(1) = I_1;
R_RK2_3(1) = R_1;

% Step 2: Loops Over All Time Points For h1

% Step 2.1: Loop For NSFDM

for j = 1:1:length(t1)-1
  S_NSFDM1(j+1) = S_NSFDM1(j)/(1+h1*alpha*I_NSFDM1(j)/N);
  I_NSFDM1(j+1) = (I_NSFDM1(j)+h1*alpha*S_NSFDM1(j+1)*I_NSFDM1(j)/N)/(1+h1*beta);
  R_NSFDM1(j+1) = R_NSFDM1(j)+h1*beta*I_NSFDM1(j+1);
endfor

% Step 2.2: Loop For EE

for j = 1:1:length(t1)-1
  S_EE1(j+1) = S_EE1(j)-h1*alpha*S_EE1(j)*I_EE1(j)/N;
  I_EE1(j+1) = I_EE1(j)+h1*alpha*S_EE1(j)*I_EE1(j)/N-h1*beta*I_EE1(j);
  R_EE1(j+1) = R_EE1(j)+h1*beta*I_EE1(j);
endfor

% Step 2.3: Loop For RK2

for j = 1:1:length(t1)-1
  k1 = [-alpha*I_RK2_1(j)*S_RK2_1(j)/N; alpha*I_RK2_1(j)*S_RK2_1(j)/N-beta*I_RK2_1(j); beta*I_RK2_1(j)];
  k2 = [-alpha*(I_RK2_1(j)+h1/2*alpha*I_RK2_1(j)*S_RK2_1(j)/N-h1/2*beta*I_RK2_1(j))*(S_RK2_1(j)-h1/2*alpha*I_RK2_1(j)*S_RK2_1(j)/N)/N; ...
        alpha*(I_RK2_1(j)+h1/2*alpha*I_RK2_1(j)*S_RK2_1(j)/N-h1/2*beta*I_RK2_1(j))*(S_RK2_1(j)-h1/2*alpha*I_RK2_1(j)*S_RK2_1(j)/N)/N-beta*(I_RK2_1(j)+h1/2*alpha*I_RK2_1(j)*S_RK2_1(j)/N-h1/2*beta*I_RK2_1(j)); ...
        beta*(I_RK2_1(j)+h1/2*alpha*I_RK2_1(j)*S_RK2_1(j)/N-h1/2*beta*I_RK2_1(j))];
  S_RK2_1(j+1) = S_RK2_1(j) + h1*k2(1);
  I_RK2_1(j+1) = I_RK2_1(j) + h1*k2(2);
  R_RK2_1(j+1) = R_RK2_1(j) + h1*k2(3);
endfor

% Step 3: Loops Over All Time Points For h2

% Step 3.1: Loop For NSFDM

for j = 1:1:length(t2)-1
  S_NSFDM2(j+1) = S_NSFDM2(j)/(1+h2*alpha*I_NSFDM2(j)/N);
  I_NSFDM2(j+1) = (I_NSFDM2(j)+h2*alpha*S_NSFDM2(j+1)*I_NSFDM2(j)/N)/(1+h2*beta);
  R_NSFDM2(j+1) = R_NSFDM2(j)+h2*beta*I_NSFDM2(j+1);
endfor

% Step 3.2: Loop For EE

for j = 1:1:length(t2)-1
  S_EE2(j+1) = S_EE2(j)-h2*alpha*S_EE2(j)*I_EE2(j)/N;
  I_EE2(j+1) = I_EE2(j)+h2*alpha*S_EE2(j)*I_EE2(j)/N-h2*beta*I_EE2(j);
  R_EE2(j+1) = R_EE2(j)+h2*beta*I_EE2(j);
endfor

% Step 3.3: Loop For RK2

for j = 1:1:length(t2)-1
  k1 = [-alpha*I_RK2_2(j)*S_RK2_2(j)/N; alpha*I_RK2_2(j)*S_RK2_2(j)/N-beta*I_RK2_2(j); beta*I_RK2_2(j)];
  k2 = [-alpha*(I_RK2_2(j)+h2/2*alpha*I_RK2_2(j)*S_RK2_2(j)/N-h2/2*beta*I_RK2_2(j))*(S_RK2_2(j)-h2/2*alpha*I_RK2_2(j)*S_RK2_2(j)/N)/N; ...
        alpha*(I_RK2_2(j)+h2/2*alpha*I_RK2_2(j)*S_RK2_2(j)/N-h2/2*beta*I_RK2_2(j))*(S_RK2_2(j)-h2/2*alpha*I_RK2_2(j)*S_RK2_2(j)/N)/N-beta*(I_RK2_2(j)+h2/2*alpha*I_RK2_2(j)*S_RK2_2(j)/N-h2/2*beta*I_RK2_2(j)); ...
        beta*(I_RK2_2(j)+h2/2*alpha*I_RK2_2(j)*S_RK2_2(j)/N-h2/2*beta*I_RK2_2(j))];
  S_RK2_2(j+1) = S_RK2_2(j) + h2*k2(1);
  I_RK2_2(j+1) = I_RK2_2(j) + h2*k2(2);
  R_RK2_2(j+1) = R_RK2_2(j) + h2*k2(3);
endfor

% Step 4: Loops Over All Time Points For h3

% Step 4.1: Loop For NSFDM

for j = 1:1:length(t3)-1
  S_NSFDM3(j+1) = S_NSFDM3(j)/(1+h3*alpha*I_NSFDM3(j)/N);
  I_NSFDM3(j+1) = (I_NSFDM3(j)+h3*alpha*S_NSFDM3(j+1)*I_NSFDM3(j)/N)/(1+h3*beta);
  R_NSFDM3(j+1) = R_NSFDM3(j)+h3*beta*I_NSFDM3(j+1);
endfor

% Step 4.2: Loop For EE

for j = 1:1:length(t3)-1
  S_EE3(j+1) = S_EE3(j)-h3*alpha*S_EE3(j)*I_EE3(j)/N;
  I_EE3(j+1) = I_EE3(j)+h3*alpha*S_EE3(j)*I_EE3(j)/N-h3*beta*I_EE3(j);
  R_EE3(j+1) = R_EE3(j)+h3*beta*I_EE3(j);
endfor

% Step 4.3: Loop For RK2

for j = 1:1:length(t3)-1
  k1 = [-alpha*I_RK2_3(j)*S_RK2_3(j)/N; alpha*I_RK2_3(j)*S_RK2_3(j)/N-beta*I_RK2_3(j); beta*I_RK2_3(j)];
  k2 = [-alpha*(I_RK2_3(j)+h3/2*alpha*I_RK2_3(j)*S_RK2_3(j)/N-h3/2*beta*I_RK2_3(j))*(S_RK2_3(j)-h3/2*alpha*I_RK2_3(j)*S_RK2_3(j)/N)/N; ...
        alpha*(I_RK2_3(j)+h3/2*alpha*I_RK2_3(j)*S_RK2_3(j)/N-h3/2*beta*I_RK2_3(j))*(S_RK2_3(j)-h3/2*alpha*I_RK2_3(j)*S_RK2_3(j)/N)/N-beta*(I_RK2_3(j)+h3/2*alpha*I_RK2_3(j)*S_RK2_3(j)/N-h3/2*beta*I_RK2_3(j)); ...
        beta*(I_RK2_3(j)+h3/2*alpha*I_RK2_3(j)*S_RK2_3(j)/N-h3/2*beta*I_RK2_3(j))];
  S_RK2_3(j+1) = S_RK2_3(j) + h3*k2(1);
  I_RK2_3(j+1) = I_RK2_3(j) + h3*k2(2);
  R_RK2_3(j+1) = R_RK2_3(j) + h3*k2(3);
endfor

% Step 5: Plots

figure(1)
plot(t1,S_NSFDM1)
hold on
%plot(t,S_EE)
%hold on
plot(t1,S_RK2_1)
title('Susceptible Population, Time Step Size: h = 7.5')
xlabel('t')
ylabel('S(t)')
legend('NSFDM', 'RK2')
hold off

figure(2)
plot(t2,S_NSFDM2)
hold on
plot(t2,S_RK2_2)
hold on
plot(t2,S_EE2)
title('Susceptible Population, Time Step Size: h = 5')
xlabel('t')
ylabel('S(t)')
legend('NSFDM', 'RK2', 'EE')
hold off

figure(3)
plot(t3,S_NSFDM3)
hold on
plot(t3,S_RK2_3)
hold on
plot(t3,S_EE3)
title('Susceptible Population, Time Step Size: h = 0.1')
xlabel('t')
ylabel('S(t)')
legend('NSFDM', 'RK2', 'EE')
hold off

figure(4)
plot(t3,S_NSFDM3+I_NSFDM3+R_NSFDM3)
hold on
plot(t3,S_RK2_3+I_RK2_3+R_RK2_3)
hold on
plot(t3,S_EE3+I_EE3+R_EE3)
title('Total Population Size Conservation, Time Step Size: h = 0.1')
xlabel('t')
ylabel('N')
legend('NSFDM','RK2','EE')
hold off