function[s,result,counter]=quasi_newton(cx,cy,xi,yi,epsilon_newton) % d est un doublet cx, cy
    dbtype('grad_ctls.m') ;
    dbtype('ctls.m') ;
    dbtype('fletcher.m')
    % On initialise la matrice Hk à 0 
    % Attention dans ce code je choisit la notation Hk mais ici en réalité
    % nous avons Hk^-1 c'est la matrice inverse de Hk

    I = eye(2);
    alpha = 1 ;
    Hk = eye(2); % On doit avoir une matrice de taille 2
    disp(Hk)
    result = [] ;
    counter = 0 ;


    Grad_fk = grad_ctls(cx,cy,xi,yi)

    while (norm(Grad_fk)>epsilon_newton)
        Grad_fk = grad_ctls(cx,cy,xi,yi)
        dk = -Hk*Grad_fk ;   
        
        alpha = fletcher(cx,cy,alpha,dk,xi,yi);
        
        cx = cx + alpha*dk(1) ;
        cy = cy + alpha*dk(2) ;
    
        Grad_fk1 = grad_ctls(cx,cy,xi,yi) ;
    
        yk = Grad_fk1 - Grad_fk ; 
    
        dk_ = alpha*dk ; % dk barre On peut aussi l'écrire comme la différence des xk+1 xk
    
        Hk = (I - dk_*(yk')/((dk_')*yk) )*Hk*(I- yk*(dk_')/((dk_')*yk)  ) + dk_*(dk_')/((dk_')*yk) ;

        result = [result;[cx,cy]] ;
        counter = counter + 1
    end
    s = [cx,cy] ;

end
