function[d]=grad_ctls(cx,cy,xi,yi)
    n = length(xi) ; 
    d = [0;0] ;
    R = 1.5 ;

    %syms k ;
    %d = symsum( ( sqrt((xi(k)-cx)^2 + (yi(k) - cy)^2)-R )^2   ,k,0,29   ) ;

    for k = 1:length(xi)

        d(1) = d(1) + 2*(xi(k)-cx)*(R/(sqrt((xi(k)-cx)^2 + (yi(k) - cy)^2)) - 1 );
        d(2) = d(2) + 2*(yi(k)-cy)*(R/(sqrt((xi(k)-cx)^2 + (yi(k) - cy)^2)) - 1 );

    end
end
