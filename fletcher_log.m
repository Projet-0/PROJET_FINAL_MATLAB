function[s]=fletcher_log(cx,cy,a,d,xi,yi,sigma) % d est un doublet cx, cy
    dbtype('grad_ctls_log.m') ;
    dbtype('ctls_log.m') ;

    al = 0 ; % alpha left
    ar = inf ; % alpha right
    
    lambda = 20 ;

    B1 = 10^-3 ;
    B2 = 0.5 ;
    %y = 0 ; % : Gamma On le déclare et on lui donne sa valeur

    % Nous devons calculer gamma à partir de a,B1,B2, d est moins le
    % gradient

    y = -B1*(grad_ctls_log(cx,cy,xi,yi,sigma)') * d; % Il faut que d soit un vecteur colonne 



    % Condition de Wolf 1 et Wolf 2 Il vaut mieux faire un while
    while (ctls_log(cx+a*d(1),cy+a*d(2),xi,yi,sigma) > ctls_log(cx,cy,xi,yi,sigma) - a*y ) || ((grad_ctls_log(cx+a*d(1),cy+a*d(2),xi,yi,sigma)')*d/((grad_ctls_log(cx,cy,xi,yi,sigma)')*d) > B2  ) % On calcule la fonction de cout et on vérifie 
        
        % Condition de Wolf 1 
        if(ctls_log(cx+a*d(1),cy+a*d(2),xi,yi,sigma) <= ctls_log(cx,cy,xi,yi,sigma) - a*y ) % Si la condition de Wolf 1 est respectée
            al = a;
                if(ar >= inf )
                    a = lambda*a ;
                else
                    a = (al+ar)/2 ;
                end
        else % Si elle n'est pas respectée
            ar = a ;
            a = (al + ar)/2 ;
        end
    end
    s = a ;
end
