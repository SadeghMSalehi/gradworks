function [ Y ] = entropy_min( X )
%ENTROPY - minimize X to have minimum entropy
%   it follows gradient descent direction
[r c] = size(X);
f = @(t,Y)(gradientEntropy(t,Y,r,c));
f2 = @(t,Y)(goToMean(t,Y,r,c));
Xo = reshape(X, [r*c 1]);

[t, Xo] = ode45(f, [0 30], Xo);
% for j=1:1000
%     X = reshape(Xo, [r c]);
%     Xm = repmat(mean(X), [r, 1]);
%     Xc = X - Xm;    
%     XXt = Xc*Xc'/c+eye(r,r);
%     invXXt = inv(XXt);
%     dXXtdX = invXXt * Xc;
%     Xc = Xc - 0.1 *dXXtdX;
%     X = Xc + Xm;
%     Xo = reshape(X, [1, r*c]);
% end
G = [];
for j=1:size(Xo,1)
    X = reshape(Xo(j,:), [r c]);
    d = det(X*X');
    G(end+1) = d;
end
%t
plot(t, Xo, 'o-');
figure;
plot(t, G, 'r*-');
Y = reshape(Xo(end,:), [r c]);
end

function [ dX ] = goToMean(t, Xo, r, c)
    X = reshape(Xo, [r c]);
    Xm = repmat(mean(X), [r, 1]);
    Xc = X - Xm;    
    dX = -reshape(Xc, [r*c 1]);
end

function [ dX ] = gradientEntropy(t,Xo,r,c)
    X = reshape(Xo, [r c]);
    Xm = repmat(mean(X), [r, 1]);
    Xc = X - Xm;    
    XXt = Xc*Xc'/c+eye(r,r);
    invXXt = inv(XXt);
    dX = invXXt * Xc;
    dX = reshape(dX, [r*c 1]);
    dX = -dX;
end