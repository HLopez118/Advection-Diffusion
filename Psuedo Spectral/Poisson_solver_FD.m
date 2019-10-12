function P = Poisson_solver_FD(N,w)
% clear all
% N = 4;
KT = 2*N;
L = 1;
x = linspace(-L,L,KT+1);
x = x(1:KT);
y = x;
[X,Y] = meshgrid(x,y);

h = 2*L/(KT+1);

% w = w0(KT,x,y,X,Y);
N = KT;
% B = zeros(N,N);

% B(1,1) = -4;
% B(1,2) = 1;
% for i = 2:N-1
%     B(i,i) = -4;
%     B(i,i-1) = 1;
%     B(i,i+1) = 1;
% end
% B(end,end) = -4;
% B(end,end-1) = 1;

b = ones(N,1);

B = spdiags([b -4*b b],[-1 0 1], N,N);

B(1,end) = 1;
B(end,1) = 1;

A = eye(length(B));

T1 = kron(A,B);

% T1(1,end) = 1;
% T1(end,1) = 1;
% 
% I1 = zeros(N,N);
% I1(1,3) = 1;
% for j = 3:N-2
%     I1(j,j+2) = 1;
%     I1(j,j-2) = 1;
% end
% I1(end-2,end) = 1;
% I1(end-3,end-1) = 1;
a = zeros(N,1);
I1 = spdiags([b a b],[-1 0 1],N,N);

T3 = kron(I1,A);

I = zeros(N,N);
I(1,end) = 1;
I(end,1) = 1;

T2 = kron(I,A);


T = T1+T3+T2;

w = reshape(w,KT^2,1);

P = (1/(h^2))*w\T;

P = reshape(P,KT,KT);

% mesh(P)
end