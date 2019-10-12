clc 
clear all
dt = .5;
tmax = 30;
t = 0:dt:tmax;
tlength = length(t);
N = 4;
KT = 2*N;

v = 0.0001;

L = 1;
h = L/N;
x = linspace(-L,L,KT+1);
x = x(1:KT);
y = x;

[X,Y] = meshgrid(x,y);

w = w0(KT,x,y,X,Y);

% 

% w1 = (h^2)*w;

% B = [-1 -4 1];

% I = eye((N+1));

% A = kron(I,B);
w = reshape(w',KT^2,1);

A = zeros(KT^2,KT^2);

for n = 1:KT^2
    A(n,n) = -4;
end

for n = 2:(KT^2)-1
    if mod(n,KT) == 0
        A(n,n+1) = 0;
        A(n,n-1) = 1;
    elseif mod(n-1,KT) == 0
        A(n,n+1) = 1;
        A(n,n-1) = 0;
    else
        A(n,n+1) = 1;
        A(n,n-1) = 1; 
    end
end

A(1,2) = 1;
A(end,end-1) = 1;

for n = 4:(KT^2)-4
    if mod(n-1,KT) == 0
        A(n,n+3) = 1;
        A(n,n-3) = 0;
    elseif mod(n-4,KT) == 0
        A(n,n+3) = 0;
        A(n,n-3) = 1;
    else
        A(n,n+3) = 0;
        A(n,n-3) = 0; 
    end
end
stf = 1;
A(1,4) = 1;
A(end,end-3) = 1;
for i = 1
    stf = stf + KT*(i);
    for n = 1:KT^2 - (stf)
        A(n,n+(stf-1)) = 1;
        A(n+(stf-1),n) = 1;
    end
end



B = inv(A);

sf = (h^2)*A*w;

sf = reshape(sf.',KT,KT).';

mesh(sf)

