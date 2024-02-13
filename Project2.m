clc, clear all, close all, format compact
beep off

%% Center frequency and sampling intervals 

wc = 8.54 ;
delt = .01 ;                    % Time sampling period [sec]
T = 10.28 ;                     % Time record length [sec]
N = T/delt ;                    % Number of samples in record
ts = 9.06 ;                     % Settling Time
sim_len_min = 20*ts ;           % 181.2s
m = 18 ;                        % total number of N-periods/Tperiods simulated
t = [0:m*N-1]*delt ;            % Simulation length
delw = 2*pi/T ;                 % Frequency sampling resolution [rad/sec]
wn = pi/delt ;                  % Nyquist rate [rad/sec]
ws = 2*pi/delt ;                % Sampling rate [rad/sec]  %.1 to 50
w = [0:N/2-1]*delw ;            % FFT frequency vector up to wn
w_m = [0:N*m-1]*delw ;
w_f = [0:N-1]*delw ;
f = w/(2*pi) ;                  % FFT freq in Hz
rand('state',0)

%% Plant
Za = -140.22 ;
Zq = -1.456 ;
Ma = -2.01 ;
Mq = -.35 ;
Zde = -25.9 ;
Mde = -5 ;

A = [Za Zq; Ma Mq] ;
Aa = [Za Zq; Ma Mq] ;
B = [Zde; Mde] ;
C = [0 1] ;
D = 0 ;
G = ss(A,B,C,D) ;

[num,den] = ss2tf(A,B,C,D) ;
OL_tf = tf(num,den) ;
[mag,phase] = bode(OL_tf,w) ;
Mag = reshape(mag,1,length(w)) ;
Phase = reshape(phase,1,length(w)) ;
plant = zeros(1,length(w_m)) ;
for i = 1:length(w_m)
    plant(:,i) = evalfr(OL_tf,j*w_m(i)) ;
end

% control 
w0 = 8 ;
a = 1/100 ;
M = 1.5 ;
s = tf('s') ;
W1 = (s/M+w0)/(s+w0*a) ;
W2 = (s+.1)/(s+100) ;
W3 = [];

[K,CL,Gam,INFO] = mixsyn(G,W1,W2,W3) ;
[numK,denK] = ss2tf(K.A,K.B,K.C,K.D) ;
K_tf = tf(numK,denK) ;

[Magk,Phasek] = bode(K,w) ;
magk = reshape(Magk,1,length(w)) ;
phasek = reshape(Phasek,1,length(w)) ;

%% Noise plants
%Dryeden Turbulence
b = 27*3 ; % wingspan in m
h = 3000 ; % altitude in m
V = 385/1.944 ; % speed knots to m/s
Lw = h/2 ;
Lv = h/2/(.177+.000823*h)^(1.2) ;
W20 = 15/1.944 ;      % wind speed of 15 knots in m at height of 60m
sig_w = .1*W20 ;
s = tf('s') ;


Hw = sig_w*sqrt(2*Lv/(pi*V))*(1+2*sqrt(3)*Lw/V*s)/(1+2*Lw/V*s)^2 ;
Hq = s/V/(1+4*b/(pi*V)*s)*Hw ;
v_in = impulse(Hq,t)+.01*rand(size(t))' ;      % Simulate Dryden function with psuedo random input & noise

%% Psuedo Random Input 
%psuedo random input
upRMS = 1 ;
up0 = zeros(size(t))+sqrt(2/N)*upRMS ;  

for k = 1:N/2-1
    A = 2/sqrt(N) ;             % Vary amplitude
    phi = (2*pi)*rand ;      % Random phase uniformly dist. in 0:2pi
    up0 = up0+A*cos(2*pi*k/T*t+phi) ;
end

up0_1 = up0(:,1:N) ;
up0_3 = up0(:,2*N+1:3*N) ;
up0_m = up0((m-1)*N+1:m*N) ;

% Psuedo random input
up = up0+v_in' ;                 % Add input noise from Dryden

up_1 = up(:,1:N) ;
up_3 = up(:,2*N+1:3*N) ;
up_m = up((m-1)*N+1:m*N) ;

%% Psuedo Random Ouput 
% Simulate plant output
wd = .01*rand(size(t)) ;
yp = lsim(num,den,up,t)' + wd ;      % Output with noise
yp_1 = yp(:,1:N) ;
yp_3 = yp(:,2*N+1:3*N) ;
yp_m = yp((m-1)*N+1:m*N) ;

yp0 = lsim(num,den,up0,t)' ;          % Output w/o noise
yp0_1 = yp0(:,1:N) ;
yp0_3 = yp0(:,2*N+1:3*N) ;
yp0_m = yp0((m-1)*N+1:m*N) ;

%% Discrete Time Model
tf_ol = c2d(OL_tf,delt,'zoh') ;

theta = [-tf_ol.Denominator{1,1}(3); -tf_ol.Denominator{1,1}(2); tf_ol.Numerator{1,1}(3); tf_ol.Numerator{1,1}(2)] ;
alpha = [-tf_ol.Denominator{1,1}(3); -tf_ol.Denominator{1,1}(2)] ;
beta = [tf_ol.Numerator{1,1}(3); tf_ol.Numerator{1,1}(2)] ;

% No noise
y0n = zeros(1,m*N) ;
mu0 = zeros(m*N,1) ;
phi0_t = zeros(m*N,4) ;
for k = 1:m*N-2
    y0n(:,1) = yp0(:,1) ;
    y0n(:,2) = yp0(:,2) ;
    phi0_t(k+1,:) = [y0n(:,k) y0n(:,k+1) up0(:,k) up0(:,k+1)] ;
    y0n(:,k+2) = phi0_t(k+1,:)*theta ;
end
y0 = phi0_t*theta ;            % Shift y to y(k+1)


% Phi'Theta output sectioning
y0_1 = y0(1:N,:) ;
y0_3 = y0(2*N+1:3*N,:) ;
y0_m = y0((m-1)*N+1:m*N,:) ;

% With noise
yn = zeros(1,N) ;
phi_t = zeros(m*N,4) ;
for k = 1:m*N-2
    yn(:,1) = yp(:,1) ;
    yn(:,2) = yp(:,2) ;
    phi_t(k+1,:) = [yn(:,k) yn(:,k+1) up(:,k) up(:,k+1)] ;
    yn(:,k+2) = phi_t(k+1,:)*theta ;
end
y = phi_t*theta ;            % Shift y to y(k+1)

% Phi'Theta output sectioning
y_1 = y(1:N,:) ;
y_3 = y(2*N+1:3*N,:) ;
y_m = y((m-1)*N+1:m*N,:) ;


%% Comparison of models
figure(1)
sgtitle('Comparison of Simulation and Discrete Time Output')
subplot(2,1,1)
semilogx(f,20*log10(abs(yp0_m(1:N/2))))
hold on
semilogx(f,20*log10(abs(y0_m(1:N/2))))
legend('Psuedo Random Output','Phi^t*Theta')
title('Without Noise')
xlabel('frequency, [Hz]')
ylabel('magnitude, [dB]')

subplot(2,1,2)
semilogx(f,20*log10(abs(yp_m(1:N/2))))
hold on
semilogx(f,20*log10(abs(y_m(1:N/2))))
legend('Psuedo Random Output','Phi^t*Theta')
title('With Noise')
xlabel('frequency, [Hz]')
ylabel('magnitude, [dB]')

%% Predictor
theta_est = 1.5*rand(4,1) ;                % assign random parameter vector
yeh = phi0_t*theta_est ;                   % prediction
e_eh = phi0_t*(theta-theta_est) ;          % Equation Error


theta_h = theta ;
y_eh = phi0_t*theta_h ;
theta_t = theta-theta_h ;
e_e = phi0_t*theta_t ;           % PROVE w/Plant Param

figure(2)
sgtitle('Predictor Validation')
subplot(2,1,1)
plot(1:N,e_eh((m-1)*N+1:m*N,:))
title('Predictor w/Random Theta')
ylabel('Error Magnitude')
xlabel('Iteration, k')
subplot(2,1,2)
plot(1:N,e_e((m-1)*N+1:m*N,:))
title('Predictor w/Plant Parameters')
ylabel('Error Magnitude')
xlabel('Iteration, k')


%% Gradient Descent

% Without Noise

for k = 1:m*N-2
    theta_hat0(:,1) = theta_est ;
    mu0(k,:) = abs(1/(phi0_t(k+1,:)*phi0_t(k+1,:)')-.1) ;
    e_g0(k+2,:) = phi0_t(k+1,:)*(theta-theta_hat0(:,k)) ;
    theta_tilde0(:,k) = (theta-theta_hat0(:,k)) ;
    norm_eg0(k,:) = theta_tilde0(:,k)'*theta_tilde0(:,k) ;
    theta_hat0(:,k+1) = theta_hat0(:,k)+2*mu0(k,:)*phi0_t(k+1,:)'*e_g0(k+2,:) ;
end

figure(3)
semilogx(1:N-2,abs(norm_eg0(1:N-2)'))
title('Error Norm w/o Noise')
xlabel('Interation Index, k')
ylabel('Error Norm Magnitude (abs)')

figure(4)
semilogx(1:N,abs(e_g0(1:N,:)'))
title('Prediction Error w/o Noise')
xlabel('Interation Index, k')
ylabel('Prediction Error (abs)')

figure(5)
plot(1:N,theta_hat0(:,1:N))
title('Iteration History w/o Noise')
xlabel('Interation Index, k')
ylabel('Parameter Estimate')
legend('-Alpha0','-Alpha1','Beta0','Beta1')

% With Noise
theta_hat = zeros(4,m*N) ;
yg_h = zeros(1,m*N) ;
e_g = zeros(m*N,1) ;
norm_eg = zeros(m*N,1) ;
flag = zeros(m*N,1) ;


for k = 1:m*N-2
    theta_hat(:,1) = theta_est ;
    mu(k,:) = abs(1/(phi_t(k+1,:)*phi_t(k+1,:)')-.1) ;
    e_g(k+2,:) = phi_t(k+1,:)*(theta-theta_hat(:,k)) ;
    theta_tilde(:,k) = (theta-theta_hat(:,k)) ;
    norm_eg(k,:) = theta_tilde(:,k)'*theta_tilde0(:,k) ;
    theta_hat(:,k+1) = theta_hat(:,k)+2*mu(k,:)*phi_t(k+1,:)'*e_g(k+2,:) ;
end

k = 1:N ;
figure(6)
semilogx(k,norm_eg(1:N)')
title('Error Norm w/Noise')
xlabel('Interation Index, k')
ylabel('Error Norm Magnitude (abs)')

figure(7)
semilogx(k,e_g(1:N,1)')
title('Prediction Error w/Noise')
xlabel('Interation Index, k')
ylabel('Prediction Error (abs)')

figure(8)
plot(k,theta_hat(:,1:N))
title('Iteration History w/Noise')
xlabel('Interation Index, k')
ylabel('Parameter Estimate')
legend('-Alpha0','-Alpha1','Beta0','Beta1')

% Transfer Function Comparison
num0 = [theta_hat0(4,m*N-1) theta_hat0(3,m*N-1)] ;
den0 = [1 -theta_hat0(2,m*N-1) -theta_hat0(1,m*N-1)] ;
tf_est0 = tf(num0,den0,delt) ;

num1 = [theta_hat(4,m*N-1) theta_hat(3,m*N-1)] ;
den1 = [1 -theta_hat(2,m*N-1) -theta_hat(1,m*N-1)] ;
tf_est = tf(num1,den1,delt) ;

figure(9)
bode(tf_ol)
hold on
bode(tf_est0)
bode(tf_est)
legend('Plant','Estimated Plant w/o Noise','Estimated Plant w/Noise')

%% Batch

% Batch Least Squares
theta_b = (phi_t'*phi_t)^-1*phi_t'*y ;
rnk = rank(phi_t'*phi_t) ;      % full rank ---> persistent excitation!

% W^T*W size known
ETE1 = (y-phi_t*theta_b)'*(y-phi_t*theta_b) ;
R = (phi_t'*phi_t) ;
eige = eig(R) ;

delta = 1 ;
d_thetaa = (abs(delta*ETE1/eige(1,1)))^(1/2) ;
d_thetab = (abs(delta*ETE1/eige(2,1)))^(1/2) ;
d_thetac = (abs(delta*ETE1/eige(3,1)))^(1/2) ;
d_thetad = (abs(delta*ETE1/eige(4,1)))^(1/2) ;
d_theta_e = [d_thetaa; d_thetab; d_thetac; d_thetad] ;
theta_eset1 = d_theta_e+theta_b ;
theta_eset2 = .5*d_theta_e+theta_b ;

num_eset1 = [theta_eset1(4,1) theta_eset1(3,1)] ;
den_eset1 = [1 -theta_eset1(2,1) -theta_eset1(1,1)] ;
tf_eset1 = tf(num_eset1,den_eset1,delt) ;

num_eset2 = [theta_eset2(4,1) theta_eset2(3,1)] ;
den_eset2 = [1 -theta_eset2(2,1) -theta_eset2(1,1)] ;
tf_eset2 = tf(num_eset2,den_eset2,delt) ;

figure(10)
sgtitle('Exact size of W^T*W')
bode(tf_ol)
hold on
bode(tf_eset1)
bode(tf_eset2)
legend('Plant','Extreme','1/2')


% W^T*W + 10%
delta = .1 ;                    % actually delta^2
d_theta1 = (abs(delta*ETE1/eige(1,1)))^(1/2) ;
d_theta2 = (abs(delta*ETE1/eige(2,1)))^(1/2) ;
d_theta3 = (abs(delta*ETE1/eige(3,1)))^(1/2) ;
d_theta4 = (abs(delta*ETE1/eige(4,1)))^(1/2) ;
d_theta = [d_theta1; d_theta2; d_theta3; d_theta4] ;
theta_set1 = d_theta+theta_b ;
theta_set2 = .5*d_theta+theta_b ;

num_set1 = [theta_set1(4,1) theta_set1(3,1)] ;
den_set1 = [1 -theta_set1(2,1) -theta_set1(1,1)] ;
tf_set1 = tf(num_set1,den_set1,delt) ;

num_set2 = [theta_set2(4,1) theta_set2(3,1)] ;
den_set2 = [1 -theta_set2(2,1) -theta_set2(1,1)] ;
tf_set2 = tf(num_set2,den_set2,delt) ;

figure(11)
sgtitle('Overestimated W^T*W')
bode(tf_ol)
hold on
bode(tf_set1)
bode(tf_set2)
legend('Plant','Extreme','1/2')
