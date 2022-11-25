% Projet d'optimisation 
load('measured_points (1).mat')
xi ;
yi ;

%clear; close all; clc %clear efface le workspace, close all ferme les
close all; clc 

dbtype('ctls.m') ;
dbtype('grad_ctls.m') ;
dbtype('fletcher.m');
dbtype('quasi_newton.m');

dbtype('ctls_log.m');
dbtype('grad_ctls_log.m');
dbtype('fletcher_log.m');
dbtype('fletcher_complet_log.m');
dbtype('quasi_newton_log.m');

ctls(0.2,0.2,xi,yi) ; 

cx = linspace(-1,1,100); % Valeur prise par Cx 
cy = linspace(-1,2,100); % Ensemble de valeur prise par Cy


% Intervalle -1 4 
cx2 = linspace(-1,4,200);
cy2 = linspace(-1,4,200);


epsilon = zeros(100,100) ; % Premier intervalle 

epsilon2 = zeros(200,200) ; % Second intervalle

for j = 1:100
    for k = 1: 100
        %epsilon = [epsilon, ctls( cx(j),cy(k),xi,yi ) ]; 
        epsilon(j,k) = ctls( cx(j),cy(k),xi,yi ) ;
    end

end


for j = 1:200
    for k = 1:200
        %epsilon = [epsilon, ctls( cx(j),cy(k),xi,yi ) ]; 
        epsilon2(j,k) = ctls( cx2(j),cy2(k),xi,yi ) ;
    end
end


% figure;
%     surf( cx,cy,epsilon)
% figure;
%     surf( cx2,cy2,epsilon2)


figure;
    contour( cx,cy,epsilon,100)
    xlabel ('cx')
    ylabel ('cy')
    axis equal
figure;
    contour( cx2,cy2,epsilon2,100)
    xlabel ('cx')
    ylabel ('cy')
    axis equal




% On remarque l'existence d'un maximum local, mais nous recherchons un/des
% minimums locaux/globaux

% On étudie la fonction sur deux intervalles différents, donc pour le
% graphe n2 on observe dans les valeurs de l'écart (la fonction de cout)
% De plus le nombre de points pour l'intervalle [-1,4] est identique que
% pour l'intervalle [-1 1] donc le graphe est moins précis


    %% Question 2

[min_cx,min_cy]=find(epsilon==min(epsilon(:)));
%cx(min_cx)
%disp(" LE Minimum de cy est :")
%cy(min_cy)

% On calcule un minimum de : 17.1935 pour cx = 0.4545 et cy = 1.1515 avec
% 100points


[min_cx2,min_cy2]=find(epsilon2==min(epsilon2(:)));

disp("Le Minimum de cx2 est :")
cx2(min_cx2)
disp("Le Minimum de cy2 est :")
cy2(min_cy2)

% On calcule un minimum de l'écarrt : 14.0446 pour cx2 = 2.6869 et cy =
% 1.3737 avec 100points

figure;
    plot(xi,yi,'+')
    %plot([cx(min_cx),cy(min_cy)])
    viscircles([0.4545,1.1515],1.5)
    viscircles([2.6869,1.3737],1.5)



% Le centre obtenu au second cercle semble cohérent, dû à la présence
% d'outlayer (le point tout à droite) est excentré donc à un "poiids" plus
% grand donc contribue à décaler le minimum global


% Ppour un intervalle de taille [-1,1]*[-1,2] On a besoin de
% (2/10^-4)*(3/10^-4) donc 600 Millions de valeurs

% Ppour un intervalle de taille [-1,4]*[-1,4] On a besoin de
% (5/10^-4)*(5/10^-4) donc 2.5 Milliards de valeurs

% On sait que pour chaque couple cx,cy on doit évaluer la fonction pour R
% [0.5,2.5]. On veut une précision de 10-4 près donc on doit multiplier le nom
% nombre d'évaluation de la fonction de coût que l'on a déterminé précédemment par 20 000.

% Nous allons donc chercher à trouver des algorithmes d'optimisations pour
% réduire le nombre d'évaluations de la fonciton de coût à réaliser


%% Question 3
% Cf. feuille

%% Question 4



% 
%  [FX,FY] = gradient(ctls) ;
% 
% FX(-1,-1,xi,yi)
% FY(-1,-1,xi,yi)
% 

h = 10^(-8) ;


vtiq_cx = (ctls(-3+h,-2,xi,yi)-ctls(-3,-2,xi,yi) ) / h ;% Composante selon cx
vtiq_cy = (ctls(-3,-2+h,xi,yi)-ctls(-3,-2,xi,yi) ) / h ;% Composante selon cy
vexp = grad_ctls(-3,-2,xi,yi) ;

e_cx = (vtiq_cx - vexp(1))/vtiq_cx  % Ecart relatif 
e_cy = (vtiq_cy - vexp(2))/vtiq_cy  % Ecart relatif 


% Le calcul du gradient semble cohérent avec le calcul du gradient par
% difféerence fini


%% Question 5 

% On calcule lors de l'échantillonnage régulier

epsilon3 = zeros(100,100) ; % Premier intervalle 
grad_epsilon3_cx = zeros(100,100) ; % Gradient selon cx
grad_epsilon3_cy = zeros(100,100) ; %  GRadient selon cy


epsilon3 = zeros(200,200) ; % Second intervalle

for j = 1:100
    for k = 1: 100
        epsilon3(j,k) = ctls( cx(j),cy(k),xi,yi ) ;
        %grad_epsilon3_cx(j,k) = grad_ctls(cx(j),cy(k),xi,yi)(1) ; % C'est une matrice de vecteur normalement
        %grad_epsilon3_cy(j,k) = grad_ctls(cx(j),cy(k),xi,yi)(2) ;

        G = grad_ctls(cx(j),cy(k),xi,yi); % On doit passer par une variable intermédiaire car lorsque le tableau est en lecture il ne peut pas etre en écriture au meme moment

        grad_epsilon3_cx(j,k) = G(1);
        grad_epsilon3_cy(j,k) = G(2);

    end

end

figure;
    contour( cx,cy,epsilon,100)
    xlabel ('cx')
    ylabel ('cy')
    axis equal
    hold on 
    

    quiver(cx,cy,grad_epsilon3_cx,grad_epsilon3_cy);
    axis equal;

    hold off

% Plus le gradient

%% Question 6 Méthode de Fletcher le Maréchal

% Il faut que notre fonction gradient renvoie un vecteur colonne IMPORTANT

% ON CHERCHE LA DIRECTION Dk à xk fixé .On se servira surement du gradient_ctls
% On évalue le minimum de alpha tel que f(xk + alpha*dk) 
% On calcule xk+1

% On vérifie que le gradient évaluée en xk ||Nabla(xk)|| <= v epsilon. Si ce n'est pas le cas on recommence

cxk = 2 ; % Première approximation du minimum en cx
cyk = 2 ; % Première approximation du minimum en cy 

%a = 50 ;
b = 50 ;
d = -grad_ctls(0.5,0.5,xi,yi) 
y = 0.5
%dbtype('fletcher.m');

% Il faut tester avec des valeurs opposées au gradient
d  = -grad_ctls(0.5,0.5,xi,yi)


fletcher(0.5,0.5,10^-6,-grad_ctls(0.5,0.5,xi,yi),xi,yi)  % Calcul de alpha
%fletcher_complet(0,3,10^-3,xi,yi,10^-4)


% ON utilise FLETCHER_Le_Marechal avec les coordonnées de départ (0.5,0.5)
[b, it, i] = fletcher_complet(-0.5,-0.5,10^-3,xi,yi,10^-4)

% Question 6.1 On représente l'ens des points
figure ;
    contour(cx2,cy2,epsilon2', 100)
    hold on
    axis equal
    
    plot(it(:,1), it(:,2))
    hold off


%% Attention important j'ai epsilon 2 qui est inversé donc probablement epsilon 1 aussi


 % Rapport sous forme d'un pdf, on peut le générer avec un rapport
 % automatique matlab. Il faut pouvoir interpréter les résultats obtenus

 %% Question 6.2

% cx_itere= linspace(-1,4,1049) ;
% cy_itere= linspace(-1,4,1049) ;
%  figure;
%    contour(cx_itere,cy_itere,[it(:,1), it(:,2)]);




%  figure;
%    plot(ctls(it(:,1), it(:,2),xi,yi),i)
%

result_q6_2 = [] ;
norm_q6_3 = [] ;
dist_q6_4 = [] ;

 for (variable = 1:size(it)) 
    result_q6_2 = [result_q6_2;ctls(it(variable,1),it(variable,2),xi,yi) ]  ;
    norm_q6_3 = [norm_q6_3;norm(grad_ctls(it(variable,1),it(variable,2),xi,yi))] ; 
 end

 for(variable = 1:size(it)-1)
    dist_q6_4 = [dist_q6_4;norm([it(variable+1,1)-it(variable,1);it(variable+1,2)-it(variable,2)])] ;
 end

figure;
    plot(result_q6_2) % calcul de la fonction de cout

    hold on 
    plot(norm_q6_3) % Calcul du gradient de la fonction de cout
    
    hold on
    plot(dist_q6_4) % Calcul de la distance entre cxk,cyk et cxk+1,cyk+1

    hold off
% Le resultat obtenu avec la minimisation de la fonction de cout a l'aide
% de la methode des plus fortes pentes est cohérent avec le resultat
% calcule au debut du projet. On obtient un minimum de 17.9

% Le calcul de la norme confirme la convergence du couple (cx cy) obtenu, vers un
% minimum 

% On remarque que en prenant un Beta2 plus grand (pour l'algorithme de Fletcher le Maréchal), 
% On a moins d'itérations (car on a un alpha plus grand donc un pas plus
% grand entre chaque iterations) Cela signifie que la distance entre les
% iterations est plus grande et cela s'observe sur le Graphe (on a
% egalement une allure zigzag attendue)


%     hold on
%     plot(abs(it(:,1)-it(:,2)))
% 
%     hold off

%% Question 7

[b_7, it_7, i_7] = fletcher_complet(2,3,10^-3,xi,yi,10^-4) ; % Deuxieme test avec des conditions initiales differentes

result_q7_2 = [] ;
norm_q7_3 = [] ;
dist_q7_4 = [] ;

 for (variable = 1:size(it_7)) 
    result_q7_2 = [result_q7_2;ctls(it_7(variable,1),it_7(variable,2),xi,yi) ]  ;
    norm_q7_3 = [norm_q7_3;norm(grad_ctls(it_7(variable,1),it_7(variable,2),xi,yi))] ; 
 end

 for(variable = 1:size(it_7)-1)
    dist_q7_4 = [dist_q7_4;norm([it_7(variable+1,1)-it_7(variable,1);it_7(variable+1,2)-it_7(variable,2)])] ;
 end

figure ;
    subplot(2,1,1)
    contour(cx2,cy2,epsilon2', 100)
    hold on
    axis equal
    
    plot(it_7(:,1), it_7(:,2))
    hold off 


%figure;
    subplot(2,1,2)
    plot(result_q7_2) % calcul de la fonction de cout

    hold on 
    plot(norm_q7_3) % Calcul du gradient de la fonction de cout
    
    hold on
    plot(dist_q7_4) % Calcul de la distance entre cxk,cyk et cxk+1,cyk+1

    hold off


% En changeant les conditions initiales (le point de depart) on atteint le
% deuxieme minimum local, cela montre que l on a d autres minimums local sur
% la fonction. Ce deuxieme minimum est du a la presence d outlayer dans les
% mesures (points extreme) qui peuvent etre du a des causes externes (des erreurs de mesure)
% ce qui montre les limites de la fonction de cout car les points tres
% extremes ont un poids plus grand dans l etude de la fonction et
% principalement dans l etude du minimum



%% Question 8



[b_newton, it_newton, i_newton] = quasi_newton(0.5,0.5,xi,yi,10^-3);

[b_newton2, it_newton2, i_newton2] = quasi_newton(2.67,1.37,xi,yi,10^-3);


figure ;
    contour(cx2,cy2,epsilon2', 100)
    hold on
    axis equal
    
    plot(it_newton(:,1), it_newton(:,2))
    hold off

    % Test affichage des solutions de la fonction de ciyt
figure;
    plot(xi,yi,'+')
    %plot([cx(min_cx),cy(min_cy)])
    viscircles(b_newton,1.5)
    viscircles(b_newton2,1.5)%


%% Question 9
% Nouvelle fonction de coût
ctls_log(0.5,0.5,xi,yi,10^-3)
Grad_log = grad_ctls_log(0.5,0.5,xi,yi,10^-3) ;

%Affichons les courbes associées

Matrix_ctls_log = zeros(200,200) ; % Matrice contenant les valeurs de la fonction de cout avec sigma = 10^-3

Matrix_ctls_log_sigma1 = zeros(200,200) ; % Matrice contenant les valeurs de la fonction de cout avec sigma = 10^-3
Matrix_ctls_log_sigma2 = zeros(200,200) ; % sigma = 0.1
Matrix_ctls_log_sigma3 = zeros(200,200) ; % sigma = 10


Matrix_grad_ctls_log_x = zeros(200,200) ; % Matrice gradient selon cx
Matrix_grad_ctls_log_y = zeros(200,200) ; % Matrice gradient selon cy

quasi_newton_log(0.5,0.5,xi,yi,10^-3,10) ;

% Tracé de cercle des solutions via quasi newton avec des sigma différent
figure;
    subplot(3,1,1) % Sigma = 10^-3
    plot(xi,yi,'+') 
    viscircles(quasi_newton_log(0.5,0.5,xi,yi,10^-3,10^-3),1.5)
    viscircles(quasi_newton_log(3,3,xi,yi,10^-3,10^-3),1.5)
    title('Résultat quasi_newton_log pour sigma = 10^-3')
    
    axis equal

    subplot(3,1,2) % Sigma = 1
    plot(xi,yi,'+') 
    viscircles(quasi_newton_log(0.5,0.5,xi,yi,10^-3,1),1.5)
    viscircles(quasi_newton_log(3,3,xi,yi,10^-3,1),1.5)
    title('Résultat quasi_newton_log pour sigma = 1')
    axis equal

    subplot(3,1,3) % Sigma = 10
    plot(xi,yi,'+') 
    viscircles(quasi_newton_log(0.5,0.5,xi,yi,10^-3,10),1.5)
    title('Résultat quasi_newton_log pour sigma = 10')
    axis equal


    min3_sigma1 = ctls_log( cx2(j),cy2(k),xi,yi,10^-3) ;
    min3_sigma2 = ctls_log( cx2(j),cy2(k),xi,yi,0.1) ;
    min3_sigma3 = ctls_log( cx2(j),cy2(k),xi,yi,10) ;

%% Question 9 Different sigma représentation
for j = 1:200
    for k = 1:200
        Matrix_ctls_log_sigma1(j,k) = ctls_log( cx2(j),cy2(k),xi,yi,10^-3) ; %sigma1 = 10^-3
        Matrix_ctls_log_sigma2(j,k) = ctls_log( cx2(j),cy2(k),xi,yi,0.1) ; %sigma2 = 0.1
        Matrix_ctls_log_sigma3(j,k) = ctls_log( cx2(j),cy2(k),xi,yi,10) ; %sigma3 = 10


        % calcul du min pour sigma1
        if (ctls_log( cx2(j),cy2(k),xi,yi,10^-3) < min3_sigma1) % on cherche le minimum dans la matrice
             min3_sigma1 = (ctls_log( cx2(j),cy2(k),xi,yi,10^-3));
             min_x3_sigma1 = cx2(j);
             min_y3_sigma1 = cy2(k);
        end  

        %calcul du min pour sigma2
        if (ctls_log( cx2(j),cy2(k),xi,yi,0.1) < min3_sigma2)
             min3_sigma2 = (ctls_log( cx2(j),cy2(k),xi,yi,0.1));
             min_x3_sigma2 = cx2(j);
             min_y3_sigma2 = cy2(k);
        end 

         %calcul du min pour sigma3
        if (ctls_log( cx2(j),cy2(k),xi,yi,10) < min3_sigma3)
             min3_sigma3 = (ctls_log( cx2(j),cy2(k),xi,yi,10));
             min_x3_sigma3 = cx2(j);
             min_y3_sigma3 = cy2(k);
        end 
    end
end

figure;
    subplot(3,1,1) 
    mesh( cx2,cy2,Matrix_ctls_log_sigma1) 
    xlabel ('cx')
    ylabel ('cy')
    title('Sigma 10^-3')

    subplot(3,1,2) 
    mesh( cx2,cy2,Matrix_ctls_log_sigma2') 
    xlabel ('cx')
    ylabel ('cy')
    title('Sigma 1')

    subplot(3,1,3) 
    mesh( cx2,cy2,Matrix_ctls_log_sigma3') 
    xlabel ('cx')
    ylabel ('cy')
    title('Sigma 10')

    % Tracé des contours
figure;
    subplot(3,1,1) 
    contour( cx2,cy2,Matrix_ctls_log_sigma1',200) 
    xlabel ('cx')
    ylabel ('cy')
    axis equal
    title('Sigma 10^-3')

    subplot(3,1,2) 
    contour( cx2,cy2,Matrix_ctls_log_sigma2',200) 
    xlabel ('cx')
    ylabel ('cy')
    axis equal
    title('Sigma 1')

    subplot(3,1,3) 
    contour( cx2,cy2,Matrix_ctls_log_sigma3',200) 
    xlabel ('cx')
    ylabel ('cy')
    axis equal
    title('Sigma 10')

    % On retrace les cercles en prenant 4 conditions initiales 

figure;
    plot(xi,yi,'+')
    viscircles([min_x3_sigma1,min_y3_sigma1],1.5)
    viscircles([min_x3_sigma2,min_y3_sigma2],1.5)
    viscircles([min_x3_sigma3,min_y3_sigma3],1.5,'color','g')
    xlabel ('cx')
    ylabel ('cy')
    axis equal

%% Question 10 Reprise Question 4

% Comparaison valeur théorique et valeur expérimentale

%% Question 4

h = 10^(-8) ; % Le pas utilisé pour le calcul de différence finie


vtiq_cx_log = (ctls_log(-1+h,-1,xi,yi,0.1)-ctls_log(-1,-1,xi,yi,0.1) ) / h ;% Composante selon cx
vtiq_cy_log = (ctls_log(-1,-1+h,xi,yi,0.1)-ctls_log(-1,-1,xi,yi,0.1) ) / h ;% Composante selon cy
vexp_log = grad_ctls_log(-1,-1,xi,yi,0.1) ;

e_cx_log = abs((vtiq_cx_log - vexp_log(1))/vtiq_cx_log)  % Ecart relatif selon cx
e_cy_log = abs((vtiq_cy_log - vexp_log(2))/vtiq_cy_log)  % Ecart relatif selon cy


% Le calcul du gradient semble cohérent avec le calcul du gradient par
% On a un écart relatif très faible devant 1% Il en résulte que le résultat
% calculé par la fonction gradient est cohérent, et est très proche du
% calcul par différence fini



%% Question 10 Reprise Question 5
for j = 1:200
    for k = 1:200
        Matrix_ctls_log(j,k) = ctls_log( cx2(j),cy2(k),xi,yi,10) ;
        Grad_log = grad_ctls_log(cx2(j),cy2(k),xi,yi,1) ;
        Matrix_grad_ctls_log_x(j,k) = Grad_log(1);
        Matrix_grad_ctls_log_y(j,k) = Grad_log(2);

        
    end
end

% Attention il faut faire les transposées des matrices pour bien les
% afficher

figure;
    contour( cx2,cy2,Matrix_ctls_log',200)
    xlabel ('cx')
    ylabel ('cy')
    axis equal

    hold on  
   
    quiver(cx2,cy2,Matrix_grad_ctls_log_x',Matrix_grad_ctls_log_y');  
    axis equal;
    
    title('Affichage des vecteurs gradients de la fonction de cout log')
    hold off



figure;
    contour( cx2,cy2,Matrix_ctls_log',200)
    xlabel ('cx')
    ylabel ('cy')
    axis equal

fletcher_complet_log(0.5,0.5,10^-3,xi,yi,10^-3,10) % En dernier argument il y a sigma
quasi_newton_log(0.5,0.5,xi,yi,10^-3,10)



%% Question 10 Partie 6 Implémenter la méthode des plus fortes pentes


% ON utilise FLETCHER_Le_Marechal avec les coordonnées de départ (0.5,0.5)
% On testera avec un sigma égal a 0.1

[b_log_11, it_log_11, i_log_11] = fletcher_complet_log(-1,-1,10^-3,xi,yi,10^-4,0.1);
[b_log_14, it_log_14, i_log_14] = fletcher_complet_log(-1,4,10^-3,xi,yi,10^-4,0.1);
[b_log_41, it_log_41, i_log_41] = fletcher_complet_log(4,-1,10^-3,xi,yi,10^-4,0.1);
[b_log_44, it_log_44, i_log_44] = fletcher_complet_log(4,4,10^-3,xi,yi,10^-4,0.1);

% Question 6.1 On représente l'ens des points
figure ;
    contour(cx2,cy2,Matrix_ctls_log_sigma2', 200)
    hold on
    axis equal
    
    plot(it_log_11(:,1), it_log_11(:,2))
    plot(it_log_14(:,1), it_log_14(:,2))
    plot(it_log_41(:,1), it_log_41(:,2))
    plot(it_log_44(:,1), it_log_44(:,2))
    hold off

    % Question 6.2 reprise Evolution de la fonction de cout

    result_log_q6_2 = [] ;
    norm_log_q6_3 = [] ;
    dist_log_q6_4 = [] ;

    % On va partir du points moins 1 
 for (variable = 1:size(it_log_11)) 
    result_log_q6_2 = [result_log_q6_2;ctls_log(it_log_11(variable,1),it_log_11(variable,2),xi,yi,0.1) ]  ;
    norm_log_q6_3 = [norm_log_q6_3;norm(grad_ctls_log(it_log_11(variable,1),it_log_11(variable,2),xi,yi,0.1))] ; 
 end

 for(variable = 1:size(it_log_11)-1)
    dist_log_q6_4 = [dist_log_q6_4;norm([it_log_11(variable+1,1)-it_log_11(variable,1);it_log_11(variable+1,2)-it_log_11(variable,2)])] ;
 end

figure;
    plot(result_log_q6_2,'color','g') % calcul de la fonction de cout

    hold on 
    plot(norm_log_q6_3,'color','r') % Calcul du gradient de la fonction de cout
    
    hold on
    plot(dist_log_q6_4,'color','b') % Calcul de la distance entre cxk,cyk et cxk+1,cyk+1


    xlabel('Itérations')
    ylabel('fonction de cout (vert), gradient (rouge), distance(bleu) ')
    title('Evolution de la fonction de cout (vert), du gradient (rouge) , de la distance entre les itérations au cours de la recherche  ')
    hold off

    % Le resultat obtenu avec la minimisation de la fonction de cout a l'aide
% de la methode des plus fortes pentes est cohérent avec le resultat
% calcule au debut du projet. On obtient un minimum de 17.9

% Le calcul de la norme confirme la convergence du couple (cx cy) obtenu, vers un
% minimum 

% On remarque que en prenant un Beta2 plus grand (pour l'algorithme de Fletcher le Maréchal), 
% On a moins d'itérations (car on a un alpha plus grand donc un pas plus
% grand entre chaque iterations) Cela signifie que la distance entre les
% iterations est plus grande et cela s'observe sur le Graphe (on a
% egalement une allure zigzag attendue)


%     hold on
%     plot(abs(it(:,1)-it(:,2)))
% 
%     hold off




