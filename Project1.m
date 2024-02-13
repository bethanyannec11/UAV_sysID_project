clc, clear all, close all, format compact

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

[Magk,Phasek] = bode(K,w) ;
magk = reshape(Magk,1,length(w)) ;
phasek = reshape(Phasek,1,length(w)) ;


%% Psuedo Random Input 

%psuedo random input
upRMS = 1 ;
up0 = zeros(size(t))+sqrt(2/N)*upRMS ;  

for k = 1:N/2-1
    A = 2/sqrt(N) ;             % Vary amplitude
    phi = (2*pi)*rand ;      % Random phase uniformly dist. in 0:2pi
    up0 = up0+A*cos(2*pi*k/T*t+phi) ;
end

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


%% Simulate input and ouput noise & Error Bounds
% Simulate plant input with noise
up0rms = sqrt(dot(up0,up0)/(m*N)) ;
up = up0+v_in' ;                 % Add input noise from Dryden
% up = up0 + .01*rand(size(t)) ;

% Input w/noise
Up = 2/N*fft(up(:,1:N)) ;       % Take one time record length
Up(1) = Up(1)/2 ;               % Fix DC gain normalization
Up3 = 2/N*fft(up(:,2*N+1:3*N)) ;
Up3(1) = Up3(1)/2 ;
Upm = 2/N*fft(up(:,(m-1)*N+1:m*N)) ;
Upm(1) = Upm(1)/2 ;
uprms = sqrt(dot(up,up)/(m*N)) ;

% Input w/o noise
Up0 = 2/N*fft(up0(:,1:N)) ;       % Take one time record length
Up0(1) = Up0(1)/2 ;               % Fix DC gain normalization
Up03 = 2/N*fft(up0(:,2*N+1:3*N)) ;
Up03(1) = Up03(1)/2 ;
Upm0 = 2/N*fft(up0(:,(m-1)*N+1:m*N)) ;
Upm0(1) = Upm0(1)/2 ;

% Simulate plant output
wd = .01*rand(size(t)) ;
yp = lsim(num,den,up0,t)' + wd ;      % Output with noise
yp0 = lsim(num,den,up0,t)' ;                           % Output w/o noise

% Measurable Output w/noise
Yp = 2/N*fft(yp(:,1:N)) ;
Yp(1) = Yp(1)/2 ;
Yp3 = 2/N*fft(yp(:,2*N+1:3*N)) ;
Yp3(1) = Yp3(1)/2 ;
Ypm = 2/N*fft(yp(:,(m-1)*N+1:m*N)) ;
Ypm(1) = Ypm(1)/2 ;

% Measurable Output w/o noise
Yp0 = 2/N*fft(yp0(:,1:N)) ;
Yp0(1) = Yp0(1)/2 ;
Yp03 = 2/N*fft(yp0(:,2*N+1:3*N)) ;
Yp03(1) = Yp03(1)/2 ;
Ypm0 = 2/N*fft(yp0(:,(m-1)*N+1:m*N)) ;
Ypm0(1) = Ypm0(1)/2 ;

%% Error Bounds
vh = 2/N*fft(v_in(1:N,:))' ;
vh(1) = vh(1)/2 ;  
wh = 2/N*fft(wd(:,1:N)) ;
wh(1) = wh(1)/2 ; 

[Mag_P,phase_P] = bode(num,den,w_m) ;
mag_P = reshape(Mag_P,1,length(w_m)) ;
Eb = 20*log10(abs(vh./Up)+abs(wh./(plant(:,1:N).*Up))) ;


figure(1)
semilogx(w_m(:,1:N)/(2*pi),(Eb))
% axis([10E-2 10E1 0 1])
xlabel('frequency, [Hz]')
ylabel('Error Bound Magnitude, [dB]')
title('Error Bounds of Noise Signals')

%% Input and Output Simulation WITH Noise

% m=1 input w/noise
figure(2)
subplot(2,2,1)
plot(t(1:N),up(:,1:N))
xlabel('time, [sec]')
ylabel('amplitude of up')
title('Pseudo-Random Input w/Noise m=1')

subplot(2,2,3)
plot(f,20*log10(abs(Up(:,1:N/2))))
xlabel('frequency, [Hz]')
ylabel('up fft magnitude, [dB]')

% m=1 input w/o noise
subplot(2,2,2)
plot(t(1:N),up0(:,1:N))
xlabel('time, [sec]')
ylabel('amplitude of up')
title('Pseudo-Random Input w/o Noise m=1')

subplot(2,2,4)
plot(f,20*log10(abs(Up0(:,1:N/2))))
xlabel('frequency, [Hz]')
ylabel('up fft magnitude, [dB]')

% m=1 output w/noise
figure(3)
subplot(2,2,1)
plot(t(1:N),yp(:,1:N))
xlabel('time, [sec]')
ylabel('amplitude of yp')
title('Pseudo-Random Output w/Noise m=1')

subplot(2,2,3)
semilogx(f,20*log10(abs(Yp(:,1:N/2))))
xlabel('frequency, [Hz]')
ylabel('yp fft magnitude, [dB]')

% m=1 output w/o noise
subplot(2,2,2)
plot(t(1:N),yp0(:,1:N))
xlabel('time, [sec]')
ylabel('amplitude of yp')
title('Pseudo-Random Output w/o Noise m=1')

subplot(2,2,4)
semilogx(f,20*log10(abs(Yp0(:,1:N/2))))
xlabel('frequency, [Hz]')
ylabel('yp fft magnitude, [dB]')


% m=3 input w/noise
figure(4)
subplot(2,2,1)
plot(t(2*N+1:3*N),up(:,2*N+1:3*N))
xlabel('time, [sec]')
ylabel('amplitude of up')
title('Pseudo-Random Input w/Noise m=3')

subplot(2,2,3)
plot(f,20*log10(abs(Up3(:,1:N/2))))
xlabel('frequency, [Hz]')
ylabel('up fft magnitude, [dB]')

% m=3 input w/o noise
subplot(2,2,2)
plot(t(2*N+1:3*N),up0(:,2*N+1:3*N))
xlabel('time, [sec]')
ylabel('amplitude of up')
title('Pseudo-Random Input w/o Noise m=3')

subplot(2,2,4)
plot(f,20*log10(abs(Up03(:,1:N/2))))
xlabel('frequency, [Hz]')
ylabel('up fft magnitude, [dB]')

% m=3 output w/noise
figure(5)
subplot(2,2,1)
plot(t(2*N+1:3*N),yp(:,2*N+1:3*N))
xlabel('time, [sec]')
ylabel('amplitude of yp')
title('Pseudo-Random Output w/Noise m=3')

subplot(2,2,3)
semilogx(f,20*log10(abs(Yp3(:,1:N/2))))
xlabel('frequency, [Hz]')
ylabel('yp fft magnitude, [dB]')

% m=3 output w/o noise
subplot(2,2,2)
plot(t(2*N+1:3*N),yp0(:,2*N+1:3*N))
xlabel('time, [sec]')
ylabel('amplitude of yp')
title('Pseudo-Random Output w/o Noise m=3')

subplot(2,2,4)
semilogx(f,20*log10(abs(Yp03(:,1:N/2))))
xlabel('frequency, [Hz]')
ylabel('yp fft magnitude, [dB]')

% m=18 input w/noise
figure(6)
subplot(2,2,1)
plot(t((m-1)*N+1:m*N),up(:,(m-1)*N+1:m*N))
xlabel('time, [sec]')
ylabel('amplitude of up')
title('Pseudo-Random Input w/Noise m=18')

subplot(2,2,3)
plot(f,20*log10(abs(Upm(:,1:N/2))))
xlabel('frequency, [Hz]')
ylabel('up fft magnitude, [dB]')

% m=18 input w/o noise
subplot(2,2,2)
plot(t((m-1)*N+1:m*N),up0(:,(m-1)*N+1:m*N))
xlabel('time, [sec]')
ylabel('amplitude of up')
title('Pseudo-Random Input w/o Noise m=18')

subplot(2,2,4)
plot(f,20*log10(abs(Upm0(:,1:N/2))))
xlabel('frequency, [Hz]')
ylabel('up fft magnitude, [dB]')

% m=18 output w/noise
figure(7)
subplot(2,2,1)
plot(t((m-1)*N+1:m*N),yp(:,(m-1)*N+1:m*N))
xlabel('time, [sec]')
ylabel('amplitude of yp')
title('Pseudo-Random Output w/Noise m=18')

subplot(2,2,3)
semilogx(f,20*log10(abs(Ypm(:,1:N/2))))
xlabel('frequency, [Hz]')
ylabel('yp fft magnitude, [dB]')

% m=18 output w/o noise
subplot(2,2,2)
plot(t((m-1)*N+1:m*N),yp0(:,(m-1)*N+1:m*N))
xlabel('time, [sec]')
ylabel('amplitude of yp')
title('Pseudo-Random Output w/o Noise m=18')

subplot(2,2,4)
semilogx(f,20*log10(abs(Ypm0(:,1:N/2))))
xlabel('frequency, [Hz]')
ylabel('yp fft magnitude, [dB]')

%% EFTE of 3 different N segments
EFTE = Yp./Up ;
EFTE3 = Yp3./Up3 ;
EFTEm = Ypm./Upm ;

% test1 = 180/pi*angle(EFTE(1:N/2)) ;
% test2 = 180/pi*unwrap(angle(EFTE(1:N/2)),pi/2) ;
% 
% test3 = angle(EFTE(1:N/2)) ;
% test4 = unwrap(angle(EFTE(1:N/2))) ;
% 
% phase = zeros(length(w),1) ;
% for i = 1:length(EFTE(1:N/2))
%     phase_new(i,:) = atand(imag(EFTE(i,:))/real(EFTE(i,:))) ;
%     phase_old(i,:) = rad2deg(angle(EFTE(i,:))) ;
% end


% m=1
figure(8)
subplot(2,1,1)
semilogx(f,20*log10(abs(EFTE(1:N/2))),'k')
hold on
semilogx(f,20*log10(Mag),'k-.')
legend('ETFE','analytic')
ylabel('magnitude of EFTE, [dB]')
title('ETFE at m=1 using psuedo random input')

subplot(2,1,2)
semilogx(f,180/pi*(unwrap(angle(EFTE3(1:N/2)))),'k')
hold on
semilogx(f,Phase,'k-.')
legend('ETFE','analytic')
ylabel('angle of EFTE, [deg]')
xlabel('frequency, [Hz]')
axis([f(1) f(N/2) -180 180])

% m=3
figure(9)
subplot(2,1,1)
semilogx(f,20*log10(abs(EFTE3(1:N/2))),'k')
hold on
semilogx(f,20*log10(Mag),'k-.')
legend('ETFE','analytic')
ylabel('magnitude of EFTE, [dB]')
title('ETFE at m=3 using psuedo random input')

subplot(2,1,2)
semilogx(f,180/pi*(unwrap(angle(EFTE3(1:N/2)))),'k')
hold on
semilogx(f,Phase,'k-.')
legend('ETFE','analytic')
ylabel('angle of EFTE, [deg]')
xlabel('frequency, [Hz]')
axis([f(1) f(N/2) -180 180])


% m=18
figure(10)
subplot(2,1,1)
semilogx(w,20*log10(abs(EFTEm(1:N/2))),'k')
hold on
semilogx(w,20*log10(Mag),'k-.')
legend('ETFE','analytic')
ylabel('magnitude of EFTE, [dB]')
title('ETFE at m = 18 using psuedo random input')

subplot(2,1,2)
semilogx(w,180/pi*(unwrap(angle(EFTE3(1:N/2)))),'k')
hold on
semilogx(w,Phase,'k-.')
legend('ETFE','analytic')
ylabel('angle of EFTE, [deg]')
xlabel('frequency, [Hz]')
axis([f(1) f(N/2) -300 300])


%% Error vs. Accuracy Standard m=18

E1 = (abs(EFTE(1:N/2))-abs(EFTEm(1:N/2)))./abs(EFTEm(1:N/2))*100 ;
E3 = (abs(EFTE3(1:N/2))-abs(EFTEm(1:N/2)))./abs(EFTEm(1:N/2))*100 ;


figure(11)
semilogx(f,E1)
ylabel('EFTE error magnitude, [percent]')
xlabel('frequency, [Hz]')
axis([f(1) f(N/2) -100 100])
title('Error in ETFE Using Pseudo Random Input m=1')

figure(12)
semilogx(f,E3)
ylabel('EFTE error magnitude, [percent]')
xlabel('frequency, [Hz]')
axis([f(1) f(N/2) -50 50])
title('Error in ETFE Using Pseudo Random Input m=3')

%% EFTE Average Comparison
Up_avg = 1/8*(2/N*fft(up(:,10*N+1:11*N))+2/N*fft(up(:,11*N+1:12*N))+2/N*fft(up(:,12*N+1:13*N))+2/N*fft(up(:,13*N+1:14*N))+2/N*fft(up(:,14*N+1:15*N))+2/N*fft(up(:,15*N+1:16*N))+2/N*fft(up(:,16*N+1:17*N))+2/N*fft(up(:,17*N+1:18*N))) ;
Yp_avg = 1/8*(2/N*fft(yp(:,10*N+1:11*N))+2/N*fft(yp(:,11*N+1:12*N))+2/N*fft(yp(:,12*N+1:13*N))+2/N*fft(yp(:,13*N+1:14*N))+2/N*fft(yp(:,14*N+1:15*N))+2/N*fft(yp(:,15*N+1:16*N))+2/N*fft(yp(:,16*N+1:17*N))+2/N*fft(yp(:,17*N+1:18*N))) ;
EFTE_avg = Yp_avg./Up_avg ;

figure(13)
subplot(2,1,1)
semilogx(f,20*log10(abs(EFTE_avg(1:N/2))),'k')
hold on
% semilogx(f,20*log10(abs(EFTE(1:N/2))),'b')
% semilogx(f,20*log10(abs(EFTE3(1:N/2))),'r')
% semilogx(f,20*log10(abs(EFTEm(1:N/2))),'c')
% legend('ETFE average','EFTE m=1','EFTE m=3','EFTE m=18')
ylabel('magnitude of EFTE, [dB]')
title('ETFE average using psuedo random input')

subplot(2,1,2)
semilogx(f,180/pi*angle(EFTE_avg(1:N/2)),'k')
hold on
% semilogx(f,180/pi*angle(EFTE(1:N/2)),'b')
% semilogx(f,180/pi*angle(EFTE3(1:N/2)),'r')
% semilogx(f,180/pi*angle(EFTEm(1:N/2)),'c')
% legend('ETFE average','EFTE m=1','EFTE m=3','EFTE m=18')
ylabel('angle of EFTE, [deg]')
xlabel('frequency, [Hz]')
axis([f(1) f(N/2) -300 300])

%% EFTE Nyquist and Bode Analysis
mag_loop = abs(EFTE_avg(1:N/2)).*magk ;
phase_loop = (angle(EFTE_avg(1:N/2))+pi/180*phasek) ;
Re_loop = mag_loop.*cos(phase_loop) ;
Im_loop = mag_loop.*sin(phase_loop) ;

figure(14)
plot(Re_loop,Im_loop,'b')
hold on
plot(Re_loop,-Im_loop,'b')
nyquist(OL_tf*K)
title('Nyquist Plot of EFTE')
xlabel('Re(s)')
ylabel('Im(s)')
legend('EFTE','','Plant')

figure(15)
subplot(2,1,1)
semilogx(f,20*log10(mag_loop),'k')
ylabel('magnitude of EFTE, [dB]')
title('ETFE average using psuedo random input')

subplot(2,1,2)
semilogx(f,180/pi*phase_loop,'k')
ylabel('angle of EFTE, [deg]')
xlabel('frequency, [Hz]')
% axis([f(1) f(N/2) -300 300])

%% Parameter Fitting
U = Up_avg(1:N/2) ;
Y = Yp_avg(1:N/2) ;

n = 2 ;         % # of poles
m = 1 ;         % # of zeros

% Frequency Weighting Function
W = ones(size(w)) ;       % Equal weighting
XT = W.*(j*w).^(-n).*Y ;
Ybar = W.*Y ;
for k = 2:n
        XT = [XT; W.*(j*w).^(-n+k-1).*Y];
end
for l = 1:m+1
        XT = [XT; W.*(j*w).^(-n+l-1).*U];
end

theta = inv(real(XT(:,2:length(w))*XT(:,2:length(w))'))*(real(conj(Ybar(:,2:length(w))*XT(:,2:length(w))')))' ;

% Form transfer function
denid = [1];
numid = [];
for k = 1:n
    denid = [denid, -theta(n-k+1)];
end
for l = 1:m+1
    numid = [numid, theta(n+m-l+2)];
end

ID_tf = zpk(tf(numid,denid)) ;
ID_tf = tf(numid(1),[1 .329]) ;
[MagID,PhaseID] = bode(ID_tf,w) ;
magID = reshape(MagID,1,length(w)) ;
phaseID = reshape(PhaseID,1,length(w)) ;

figure(16)
subplot(2,1,1)
semilogx(f,20*log10(abs(EFTE_avg(1:N/2))))
hold on
semilogx(f,20*log10(magID))
semilogx(f,20*log10(Mag)) 
legend('EFTE','Parametric','Analytic')
ylabel('magnitude of EFTE, [dB]')
title('Comparison of Parametric Fit to Analytic and EFTE')

subplot(2,1,2)
semilogx(f,180/pi*unwrap(angle(EFTE_avg(1:N/2))))
hold on
semilogx(f,phaseID)
semilogx(f,Phase)
legend('EFTE','Parametric','Analytic')
ylabel('angle of EFTE, [deg]')
xlabel('frequency, [Hz]')
axis([10E-2 10E1 -300 300])

%% Parametric Nyquist and Bode

figure(17)
nyquist(ID_tf*K)

figure(18)
bode(ID_tf*K)

%% Closed Loop Behavior

figure(19)
step(feedback(ID_tf*K,1))

%% Error Bound Comparison
e3 = 20*log10(abs(EFTE3(1:N))-abs(EFTEm(1:N))) ;
b = 63*ones(1,N) ;
figure(20)
plot(w_m(1:N)/(2*pi),e3)
hold on
plot(w_m(1:N)/(2*pi),20*log10(mag_P(1:N))'+b,'r') 
plot(w_m(1:N)/(2*pi),20*log10(mag_P(1:N))'-b,'r')
xlabel('Frequency, [Hz]')
ylabel('Magnitude [Db]')
legend('EFTE Error', 'P+deltaP','')

%% sigma and rho
delta = .5 ;
epsilon = .05 ;

[p1P0Q,phase1] = bode(1+OL_tf*K,w) ;
magp1P0Q = reshape(p1P0Q,1,length(w)) ;

[Q,phaseq] = bode(K,w) ;
magQ = reshape(Q,1,length(w)) ;

[PQ,phasepq] = bode(OL_tf*K,w) ;
magPQ = reshape(PQ,1,length(w)) ;

sig_analytic = 20*log10((magp1P0Q-delta)./magQ) ;

rho_analytic = 20*log10((magPQ-1/epsilon)./magQ) ;

% Parametric
[p1PQ_param,phaseparam] = bode(1+ID_tf*K,w) ;
magp1PQ_param = reshape(p1PQ_param,1,length(w)) ;

[PQ_param,phasepq_param] = bode(ID_tf*K,w) ;
magPQ_param = reshape(PQ_param,1,length(w)) ;

sig_param = 20*log10((magp1PQ_param-delta)./magQ) ;

rho_param = 20*log10((magPQ_param-1/epsilon)./magQ) ;

figure(21)
semilogx(w,sig_analytic,w,sig_param,'--')
title('Delta-Stability Function')
xlabel('frequency [Hz]')
ylabel('Magnitude [dB]')
legend('Analytic','Parametric')


figure(22)
semilogx(w,rho_analytic,w,rho_param,'--')
title('Epsilon-Performance Function')
xlabel('frequency [Hz]')
legend('Analytic','Parametric')
ylabel('Magnitude [dB]')