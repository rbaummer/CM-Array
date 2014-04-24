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
h = scatterplot(w_init'*x,1,0,'r.');
hold on;
scatterplot(w_end'*x,1,0,'b.',h)
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
h = scatterplot(w_init'*x,1,0,'r.');
hold on;
scatterplot(w_end'*x,1,0,'b.',h)
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

