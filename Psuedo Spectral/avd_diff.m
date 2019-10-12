function dwdtF = avd_diff(t,w,dummy,v,K,K2,RK,KX,KY)

wt = reshape(w,RK,RK);
psit = -wt./K;

psix = real(ifft2(1i*KX.*psit));
psiy = real(ifft2(1i*KY.*psit));
wy = real(ifft2(1i*KY.*wt));
wx = real(ifft2(1i*KX.*wt));

% dwdt = v*(Dx2+Dy2).*w - (((Dx.*sf).*(Dy.*w))-((Dy.*sf).*(Dx.*w)));

dwdtF = -v*K2.*w+reshape(fft2(wx.*psiy-wy.*psix),RK^2,1);

% dwdtF = reshape(fft2(dwdt)',length(Dx),1);

end