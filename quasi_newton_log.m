function[s,result,counter]=quasi_newton_log(cx,cy,xi,yi,epsilon_newton,sigma) % d est un doublet cx, cy
    dbtype('grad_ctls_log.m') ;
    dbtype('ctls_log.m') ;
    dbtype('fletcher_log.m') ; 

    I = eye(2);
    alpha = 1 ;
    Hk = eye(2); % On doit avoir une matrice de taille 2
    disp(Hk)
    result = [] ;
    counter = 0 ;


    Grad_fk = grad_ctls_log(cx,cy,xi,yi,sigma) ;

    while (norm(Grad_fk)>epsilon_newton)
        Grad_fk = grad_ctls_log(cx,cy,xi,yi,sigma) ;
        dk = -Hk*Grad_fk ;   
        
        alpha = fletcher_log(cx,cy,alpha,dk,xi,yi,sigma);
        
        cx = cx + alpha*dk(1) ;
        cy = cy + alpha*dk(2) ;
    
        Grad_fk1 = grad_ctls_log(cx,cy,xi,yi,sigma) ;
        yk = Grad_fk1 - Grad_fk ; 
        dk_ = alpha*dk ; 
        Hk = (I - dk_*(yk')/((dk_')*yk) )*Hk*(I- yk*(dk_')/((dk_')*yk)  ) + dk_*(dk_')/((dk_')*yk) ;
        result = [result;[cx,cy]] ;
        counter = counter + 1 ;
    end
    s = [cx,cy] ;

end
