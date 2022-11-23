

%function[s]=fletcher_function(cx,cy,xi,yi,epsilon,cxo,cyo)
function[s]=fletcher(cx,cy,a,d,xi,yi) % d est un doublet cx, cy
    dbtype('grad_ctls.m') ;
    dbtype('ctls.m') ;

    al = 0 ; % alpha left
    ar = inf ; % alpha right
    
    lambda = 20 ;

    B1 = 10^-3 ;
    B2 = 0.5 ;
    %y = 0 ; % : Gamma On le déclare et on lui donne sa valeur

    % Nous devons calculer gamma à partir de a,B1,B2, d est moins le
    % gradient

    y = -B1*(grad_ctls(cx,cy,xi,yi)') * d; % Il faut que d soit un vecteur colonne 



    % Condition de Wolf 1 et Wolf 2 Il vaut mieux faire un while
    while (ctls(cx+a*d(1),cy+a*d(2),xi,yi) > ctls(cx,cy,xi,yi) - a*y ) || ((grad_ctls(cx+a*d(1),cy+a*d(2),xi,yi)')*d/((grad_ctls(cx,cy,xi,yi)')*d) > B2  ) % On calcule la fonction de cout et on vérifie 
        
        % Condition de Wolf 1 
        if(ctls(cx+a*d(1),cy+a*d(2),xi,yi) <= ctls(cx,cy,xi,yi) - a*y ) % Si la condition de Wolf 1 est respectée
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
