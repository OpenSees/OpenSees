function processMotion()
% reshapes ground motion and creates velocity time history for analysis

% load NGA database file
rawMotion = load('GilroyNo1EW.out');
[m,n] = size(rawMotion);
nStep = m*n;

% reshape into a column vector
a = zeros(nStep,1);
k = 1;
for i = 1:m
    for j = 1:n
        a(k,1) = rawMotion(i,j);
        k = k + 1;
    end
end

% time step specified for ground motion
dT = 0.005;

% change units to m/s^2
aMod = a*9.81;

% integrate to obtain velocity time history
v = cumtrapz(aMod);
v = dT*v;

save -ascii velocityHistory.out v

return
