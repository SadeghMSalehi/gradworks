%%
% y0: initial value
% t: time steps
% f: y' = f(t,y)
function [tx,Kx] = rkf45(y0, t0, dt, t1)

% number of equations
% 3 for the homework
nVar = length(y0);

% ensure y0 is 1x3 row vector
y0 = reshape(y0, 1, nVar);

% set up for function f(t,y)
A = diag([3,1,0.01]);
f = @(t,y)(1+y.*y);

% initial value
y = y0;

% time step iteration
% currently fixed time step

% there is no spatial step
% h is equal to time step
h = dt;

% limit hMin and hMax to avoid too small or large time stpes
hMin = h/64;
hMax = 64*h;


% coefficient for h
c = [ 1/4 3/8 12/13 1 1/2 ];

% coefficient for y,k1,k2,...,k6
% kn is increment by integrating from y(t0) to y(t1)
g = [
    1/4 0 0 0 0; 
    3/32 9/32 0 0 0; 
    1932/2197 -7200/2197 7296/2197 0 0;
    439/216 -8 3680/513 -845/4104 0;
    -8/27 2 -3544/2565 1859/4104 -11/40 ];

% coefficients for RK4
yCoeff = [ 25/216 0 1408/2565 2197/4104 -1/5 0 ];

% coefficients for RK5
zCoeff = [ 16/135 0 6656/12825 28561/56430 -9/50 2/55 ];

% error coefficients from RK5 and RK4
eCoeff = zCoeff - yCoeff;

% epsilon for error
eps = 1e-5;

% iterate over from t0 to t1
% y is also function of t and t0, t1 also can be x range
t = t0;


% K is a 6xN matrix
K = zeros(6,nVar);
tx = [];
yx = [];
err = 0;
s = 1;
tx = [t err s h y y ];

% accumulation of K and error for debug
Kx = zeros(1,6);

while (t < t1)
    % next time step
    t = t + h;

    K(1,:) = h*f(t,y);
    for j = 2:6
        K(j,:) = h*f(t+c(j-1)*h,y+g(j-1,:)*K(1:5,:));
    end
    Kx(end+1,:) = K(:,1)';
    
    % y is sum of y and linear combination of Kn with yCoeff
    % K is 6x3 matrix and yCoeff is 1x6 vector
    % yCoeff*K' is 1x3 vector
    yn = y + yCoeff*K;
    
    % error between RK4 and RK5
    err = norm(eCoeff*K);
    
    % find scalar for adaptive step size
    if err == 0
        s = 0;
    else
        s = ((h/2*1e-8)/(err))^(1/4);
    end
    

    % prepare for next step
    % the large the sMin, more senstively control time step, which will
    % increase the number of iterations
    sMin = 0.75;
    
    % the smaller the sMax, more permively increase time step, which will
    % decrease the number of iterations
    sMax = 1.75;

    % 
    % if scalar is large compared to tolerance multiply half
    % if scalar is smaller, multiply two
    if (s < sMin && .5*h > hMin)
        h = .5*h;
    elseif (s > sMax && 2*h < hMax)
        h = 2*h;
    end
    
    % accumulation of t, y, err
    tx(end+1,:) = [t err s h y yn ];
        
    % y_n+1 = yn
    y = yn;
end

end