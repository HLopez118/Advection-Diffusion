function w = w0(N,x,y,xx,yy)

xn = length(x);
yn = length(y);

mpx = x(round(xn/2));
mpy = y(round(yn/2));
% h = 0.025;
yl = 10;
xl = 40*yl; %changes width of gaussian hump

a = 1;


w = a*exp(-((xl*(xx).^2)+(yl*(yy).^2)));
% d = .2;
% w = a*exp(-((xl*(xx+mpx-d).^2)+(yl*(yy+mpy-d).^2)))+a*exp(-((xl*(xx+mpx+d).^2)+(yl*(yy+mpy+d).^2)));


end
