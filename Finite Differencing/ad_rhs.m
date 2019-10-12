function dwdt = ad_rhs(t,w,dummy,h,v,T)
    
    dwdt = (v/(h^2))*T*w;
end