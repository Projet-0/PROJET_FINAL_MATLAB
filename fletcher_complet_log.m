function[s,result,counter]=fletcher_complet_log(cx,cy,a,xi,yi,epsilon_fletcher,sigma) % d est un doublet cx, cy
    dbtype('fletcher_log.m') ;
    dbtype('grad_ctls_log.m')

    s_a = [a];
    counter = 0;
    result = [cx,cy];
    % On doit d'abord calculer Dk
    while (norm(grad_ctls_log(cx,cy,xi,yi,sigma)) > epsilon_fletcher)

        d = -grad_ctls_log(cx,cy,xi,yi,sigma) ; % Page 13/30
        a = fletcher_log(cx,cy,a,d,xi,yi,sigma) ;

        
        cx = cx + a*d(1) ;
        cy = cy + a*d(2) ;

        s_a = [s_a,a] ; % Les valeurs des différents alpha
        result = [result;[cx,cy]] ;


        counter = counter + 1 ;

    end
    
    s = [cx,cy] ;
    
end
