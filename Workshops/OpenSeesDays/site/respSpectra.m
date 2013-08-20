function [p, umax, vmax, amax] = respSpectra(a, T, nStep)
% this function builds response spectra from acceleration time history

% add initial zero value to acceleration and change units
a = [0;a]*9.81;

% number of periods at which spectral values are to be computed
nPeriod = 100;

% define range of considered periods by power of 10
minPower = -3;
maxPower = 1;
% create vector of considered periods
p = logspace(minPower,maxPower,nPeriod);

% incremental circular frequency
dw = 2*pi/T;
% vector of circular freq
w = 0:dw:(nStep*dw);

% fast fourier Horm of acceleration
Afft = fft(a);

% arbitrary stiffness value
k = 1000;
% damping ratio
damp = 0.05;

% loop to compute spectral values at each period
for j = 1:nPeriod

  % compute mass and dashpot coeff to produce desired periods
    m = ((p(j)/(2*pi))^2)*k;
    c = 2*damp*(k*m)^0.5;

  % compute transfer function 
    for l = 1:(nStep/2+1)
        H(l) = 1./(-m*w(l)*w(l) + i*c*w(l) + k);
      % mirror image of Her function
        H(nStep+3-l) = conj(H(l));
    end

  % compute displacement in frequency domain using Her function
    Qfft = -m*Afft;
    for l = 1:(nStep+1)
         u(l) = H(l)*Qfft(l);
    end

  % compute displacement in time domain (ignore imaginary part)
    utime = real(ifft(u));

   % spectral displacement, velocity, and acceleration
    umax(j) = max(abs(utime));
    vmax(j) = (2*pi/p(j))*umax(j);
    amax(j) = (2*pi/p(j))*vmax(j)/9.81;
end

return
