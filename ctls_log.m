function[cout] = ctls_log(cx,cy,xi,yi,sigma)
    cout = 0 ;
    R = 1.5 ;
    for k = 1:30 
        D = sqrt((xi(k)-cx)^2 + (yi(k) - cy)^2) ; % Calcul de Di
        cout = cout + 0.5*(log(1+((D - R)^2)/(sigma^2)));
    end
end