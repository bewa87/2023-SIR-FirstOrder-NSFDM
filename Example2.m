%
% Example 2: COVID-19 Data From Spain (01 March 2020 to 15 May 2020)
%

% Step 0: Preparations

N   = 47370000;

fid = fopen('Spain-March2020-May2020.csv', 'r');
C   = textscan(fid,'%f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32', 'delimiter', ';');
fclose(fid);

spain_conf   = zeros(length(C),1);
spain_dead   = zeros(length(C),1);
spain_recov  = zeros(length(C),1);

for j = 1:1:length(C)
  spain_conf(j)  = C{1,j}(1);
  spain_dead(j)  = C{1,j}(2);
  spain_recov(j) = C{1,j}(3);
endfor

S = N - (spain_conf + spain_dead + spain_recov);
I = spain_conf - (spain_dead + spain_recov);
R = spain_dead + spain_recov;

alpha = zeros(length(C),1);
beta  = zeros(length(C),1);
t     = 1:1:length(C);
t2    = 1:0.1:length(C);

% Step 1: Estimation Of Alpha and Beta

for j = 1:1:length(C)-1
  alpha(j+1) = N/(S(j+1)*I(j))*(S(j)-S(j+1));
  beta(j+1)  = 1/(I(j+1))*(R(j+1)-R(j));
endfor

figure(1)
plot(t,I,'+')
hold on
plot(t,R,'+')
hold off
title('Infected And Recovered Population Sizes')
xlabel('t')
ylabel('I(t) And R(t)')
legend('Data I(t)','Data R(t)')

figure(2)
plot(t,alpha,'+')
hold on
plot(t2,0.52*exp(-0.032*t2),'-')
hold off
title('Time-Varying Transmission Rate')
xlabel('t')
ylabel('alpha')
legend('Data','Model: 0.52*exp(-0.03*t)')

figure(3)
plot(t,beta,'+')
title('Time-Varying Recovery Rate')
xlabel('t')
ylabel('beta')
legend('Data')

alpha_func = @(t) 0.52*exp(-0.032*t);
beta_func  = @(t) 0.045;

% Step 2: Short-Time Simulation

h          = 0.75;
t_calc     = 0:h:length(t);
S_NSFDM    = zeros(length(t_calc),1);
I_NSFDM    = zeros(length(t_calc),1);
R_NSFDM    = zeros(length(t_calc),1);
S_NSFDM(1) = N - I(1) - R(1);
I_NSFDM(1) = I(1);
R_NSFDM(1) = R(1);

% Step 2.1: Loop For NSFDM

for j = 1:1:length(t_calc)-1
  S_NSFDM(j+1) = S_NSFDM(j)/(1+h*alpha_func(j+1)*I_NSFDM(j)/N);
  I_NSFDM(j+1) = (I_NSFDM(j)+h*alpha_func(j+1)*S_NSFDM(j+1)*I_NSFDM(j)/N)/(1+h*beta_func(j+1));
  R_NSFDM(j+1) = R_NSFDM(j)+h*beta_func(j+1)*I_NSFDM(j+1);
endfor

% Step 3: Plots

figure(4)
plot(t,I,'+')
hold on
plot(t,R,'+')
hold on
plot(t_calc,I_NSFDM)
hold on
plot(t_calc,R_NSFDM)
hold off
title('Infected And Recovered Population Sizes')
xlabel('t')
ylabel('I(t) And R(t)')
legend('Data I(t)', 'Data R(t)', 'Model I(t), h = 0.75','Model R(t), h = 0.75')