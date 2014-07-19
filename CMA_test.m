function CMA_test

close all

%Angles of Arrival for QPSK signals
angles = pi/180*[20 45 74 107 137];
%Signal Power
sigma_s = [1 1 1 1 1];

%Run length (Number QPSK signals * 2)
len = 10000;
%Array length
N = 10;
%Noise power
sigma_n = 0.05; %26 dB SNR

%% generate test QPSK sources arriving at N element array
[x, V_qpsk] = QPSK_test_sources(len, N, angles, sigma_s, sigma_n);

%initial signal powers
abs(1/N*ones(1,N)*V_qpsk).^2 %#ok<NOPRT>

%% test CMA algorithm with QPSK
%set intitial weights to CBF
w_init = 1/N*ones(N,1);
%w_init = [0 0 0 0 0 1 0 0 0 0]';

%Calculate step size
%mu from "Convergence Behavior of the Constant Modulus Algorithm"
R = V_qpsk*V_qpsk';
[~,D] = eig(R);
mu = 0.25/(6*real(max(max(D))^2));

%Run CMA
[w, err] = CMA(w_init, mu, x);
w_end = w(:,end);

%check to make sure algorithm converged
assert(isfinite(w_end(1)), sprintf('Error: CMA: QPSK did not converge mu = %f',mu)); 

%plot before and after scatter plot
h = scatterplot(w_init'*x(:,2:end),1,0,'r.');
hold on;
scatterplot(w_end'*x(:,2:end),1,0,'b.',h)
title('CMA Signal Constellation');
legend('Before','After');

%plot beam pattern after convergence
figure;
bf_plot(w_end, angles);
title(sprintf('QPSK CMA Beam Pattern\n WNG: %2.1f', 1/(w_end'*w_end)));
xlabel('Degrees')
ylabel('Power (dB)');
figure;
plot(abs(err));
title('CMA Learning Curve');
xlabel('Iteration');
ylabel('|Error|');

%% test linearly constrained CMA algorithm
w_init = 1/(N-1)*ones(N-1,1);

%constraint
%sensor spacing d = lambda/2
pz = ((0:N-1) - (N-1)/2) * 0.5;
%angles in kz space, lambda normalized
lambda = 1;
kz = -2*pi./lambda*cos(angles(2));
%replica vectors
C = exp(-1i*pz'*kz);

[w, w_c, B, err] = LCCMA(w_init, mu, x, C);
w_end = w(:,end);

%GSC output
y = w_c'*x(:,2:end) - w_end'*B'*x(:,2:end);

%plot before and after scatter plot
h = scatterplot(1/N*ones(1,N)*x(:,2:end),1,0,'r.');
hold on;
scatterplot(y,1,0,'b.',h)
title('LCCMA Signal Constellation');
legend('Before','After');

%plot beam pattern after convergence
figure;
w_gsc = w_c - B*w_end;
bf_plot(w_gsc, angles);
title(sprintf('QPSK LCCMA Beam Pattern\n WNG: %2.1f', 1/(w_gsc'*w_gsc)));
xlabel('Degrees')
ylabel('Power (dB)');
figure;
plot(abs(err));
title('LCCMA Learning Curve');
xlabel('Iteration');
ylabel('|Error|');


%% test orthogonal CMA algorithm with QPSK
%set initial weights to CBF
w_init = 1/N*ones(N,1);
%w_init = [0 0 0 0 0 1 0 0 0 0]';
%set step size
mu = 0.01;
%set forgetting factor
alpha = 1-0.985;
%set initial inverse correlation matrix
R = diag(ones(N,1));

%Run o-CMA
[w, err, R_inv] = oCMA(w_init, R, mu, alpha, x);
w_end = w(:,end);

%check to make sure algorithm converged
assert(isfinite(w_end(1)), 'Error: o-CMA: QPSK did not converge'); 

%plot before and after scatter plot
h = scatterplot(w_init'*x(:,2:end),1,0,'r.');
hold on;
scatterplot(w_end'*x(:,2:end),1,0,'b.',h)
title('o-CMA Signal Constellation');
legend('Before','After');

%plot beam pattern after convergence
figure;
bf_plot(w_end, angles);
title(sprintf('QPSK O-CMA Beam Pattern\n WNG: %2.1f', 1/(w_end'*w_end)));
xlabel('Degrees')
ylabel('Power (dB)');
figure;
plot(abs(err));
title('O-CMA Learning Curve');
xlabel('Iteration');
ylabel('|Error|');

%% Test Linearly Constrained orthogonalized CMA
%set initial weights to CBF
w_init = 1/(N-1)*ones(N-1,1);
%set step size
mu = 0.01;
%set forgetting factor
alpha = 1-0.985;
%set initial inverse correlation matrix
R = diag(ones(N-1,1));

%Run LCo-CMA
[w, w_c, B, err] = LCoCMA(w_init, R, mu, x, alpha, C);
w_end = w(:,end);

%GSC output
y = w_c'*x(:,2:end) - w_end'*B'*x(:,2:end);

%plot before and after scatter plot
h = scatterplot(1/N*ones(1,N)*x(:,2:end),1,0,'r.');
hold on;
scatterplot(y,1,0,'b.',h)
title('LCO-CMA Signal Constellation');
legend('Before','After');

%plot beam pattern after convergence
w_gsc = w_c - B*w_end;
figure;
bf_plot(w_gsc, angles);
title(sprintf('QPSK LCO-CMA Beam Pattern\n WNG: %2.1f', 1/(w_gsc'*w_gsc)));
xlabel('Degrees')
ylabel('Power (dB)');
figure;
plot(abs(err));
title('LCO-CMA Learning Curve');
xlabel('Iteration');
ylabel('|Error|');

