function bf_plot(w, angles, qam_angles)

nFFT = 2048;
%Standard Array
d = 0.5;
lambda = 1;

psi = 2*pi*linspace(0,1-1/nFFT,nFFT);
mid = round((nFFT)/2)+1;
%subtract 2*pi to set axis to -pi/ to pi
psi(mid:end) = psi(mid:end) - 2*pi;
psi = fftshift(psi);

%calculate the wavenumber kz
kz = -1/d*psi;
%calculate uz space
uz = -lambda/(2*pi)*kz;

bp = 10*log10(abs(fftshift(fft(w,nFFT))).^2);
plot(acos(uz)*180/pi, bp);

%plot angles of arrival if specified
if nargin > 1
    hold on;
    for i = 1:length(angles)
        x = angles(i);
        y = (-50:0.1:0);
        plot(180/pi*x,y,'r', 'LineWidth', 2);
    end
end

if nargin > 2
    for i = 1:length(qam_angles)
        x = qam_angles(i);
        y = (-50:0.1:0);
        plot(180/pi*x,y,'c', 'LineWidth', 2);
    end
end