clc 
clear all
dt = .1;
tmax = 10;
t = 0:dt:tmax;
tlength = length(t);
N = 64;
KT = 2*N;

v = .001;

L = 1;
x = linspace(-L,L,KT+1);
x = x(1:KT);
y = x;
[X,Y] = meshgrid(x,y);

h = 2*L/(KT+1);

w = w0(KT,x,y,X,Y);

sf = poisson_sol(w,KT);

h = 2*L/(KT+1);

b = ones(KT,1);

B = spdiags([b -4*b b],[-1 0 1], KT,KT);

B(1,end) = 1;
B(end,1) = 1;

A = eye(length(B));

T1 = kron(A,B);

a = zeros(KT,1);
I1 = spdiags([b a b],[-1 0 1],KT,KT);

T3 = kron(I1,A);

I = zeros(KT,KT);
I(1,end) = 1;
I(end,1) = 1;

T2 = kron(I,A);

T = T1+T3+T2;


streamfunction = zeros(KT,KT,tlength);
streamfunction(:,:,1) = sf;
vorticity = zeros(KT,KT,tlength);
vorticity(:,:,1) = w;
tic
time = 0;
for tt = 1:tlength-1
    time = time + dt;
    Psi = poisson_sol(w,KT); % Solve Poisson Equation by SOR
         
    for i=2:KT-1
        for j=2:KT-1
            % Leap-frog scheme
            dPsi_dx(i,j) = (Psi(i+1,j)-Psi(i-1,j))/(2*h);
            dPsi_dy(i,j) = (Psi(i,j+1)-Psi(i,j-1))/(2*h);
            dw_dx(i,j) = (w(i+1,j)-w(i-1,j))/(2*h);
            dw_dy(i,j) = (w(i,j+1)-w(i,j-1))/(2*h);
            w(i,j) = w(i,j) - 2*dt*(dPsi_dy(i,j)*dw_dx(i,j) - dPsi_dx(i,j)*dw_dy(i,j));

        end
    end
    w = reshape(w,KT^2,1);
    [t1,wF1] = ode23('ad_rhs', [t(tt) t(tt+1)], w,[],h,v,T);
    w = wF1(end,:);
%     w = reshape(w,KT^2,1);
%     k1 = (1/(h^2))*T*w;
%     k2 = (1/(h^2))*T*(w+(dt*k1/2));
%     k3 = (1/(h^2))*T*(w+(dt*k2/2));
%     k4 = (1/(h^2))*T*(w+(dt*k3));
%     w = w+(dt/6)*(k1+2*(k2+k3)+k4);
    w = reshape(w,KT,KT);
    vorticity(:,:,tt) = w;
    
    
end

timer = toc;

% v = VideoWriter('Vorticity_FDS_N64.avi');
% open(v);

for j = 1:tlength-1
    pcolor(X,Y,vorticity(:,:,j))
    colorbar;
    shading flat; colormap('jet');
    title(['Advection-Diffusion, t=' num2str(t(j))])
    M(j) = getframe;
%     writeVideo(v,M(j));
end

% close(v);

% filename = 'ad_sol_FD.mat';
% save(filename,'vorticity','timer')

