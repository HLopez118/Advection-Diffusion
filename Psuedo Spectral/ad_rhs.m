function dwdt = ad_rhs(t,w,dummy,h,KT,v,T,D)
    w = reshape(w.',KT,KT).';
    sf = poisson_sol(w,KT);
    sf = reshape(sf,KT^2,1);
    w = reshape(w,KT^2,1);
    
    dwdt = (1/(2*h))*((D*sf).*(D*w)-(D*sf).*(D*w)) + v*(1/(2*h^2))*T*w;
end