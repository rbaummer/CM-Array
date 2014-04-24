%QPSK_test_sources(len, N, angles)
%returns a N element array sampling of QPSK signals arriving at angles
%specified by 'angles'.  'len' is the legnth used to generate the random
%bit stream
function [x, V] = QPSK_test_sources(len, N, angles, sigma_s, sigma_n)

%% Signal Generation
%number of QPSK signals
num = length(angles);

%generate binary bit streams for each QPSK signal
bits = randi([0 1], len, num);

qpsk = modem.oqpskmod('PhaseOffset', 0, 'SymbolOrder', 'Gray', 'InputType', 'Bit');

%generate baseband qpsk signals
s = modulate(qpsk, bits);

Nsamp = 1;

s = s.';

%% Array
%sensor spacing d = lambda/2
pz = ((0:N-1) - (N-1)/2) * 0.5;
%angles in kz space, lambda normalized
lambda = 1;
kz = -2*pi./lambda*cos(angles);
%replica vectors
V = (ones(N,1)*sigma_s.^2).*exp(-1i*pz'*kz);

%% Output
%AWGN
x = sigma_n*(randn(N,Nsamp*len)+1i*randn(N,Nsamp*len));
%sources at each sensor
x = x + V*s;

