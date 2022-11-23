function[d]=ctls(cx,cy,xi,yi)
    n = length(xi) ; 
    d = 0 ;
    R = 1.5 ;

    %syms k ;
    %d = symsum( ( sqrt((xi(k)-cx)^2 + (yi(k) - cy)^2)-R )^2   ,k,0,29   ) ;

    for k = 1:30 
        d = d + (sqrt((xi(k)-cx)^2 + (yi(k) - cy)^2)-R )^2  ;
    end
end