function[G] = grad_ctls_log(cx,cy,xi,yi,sigma)

    GradX = 0;
    GradY = 0;
    for i = 1:30
       D=sqrt((xi(i)-cx).^2+(yi(i)-cy).^2);
       gx=2*(D-1.5)/D*(cx-xi(i));
       gy=2*(D-1.5)/D*(cy-yi(i));
       GradX=GradX+1/2* gx/(sigma^2 +(D-1.5)^2);
       GradY=GradY+ 1/2* gy/(sigma^2 +(D-1.5)^2);
    end
    G = [GradX, GradY].';

end

% 
% function[G] = grad_ctls_log(cx,cy,xi,yi,sigma)
% 
%     GradX = 0;
%     GradY = 0;
%     for i = 1:30
%        D=sqrt((x(i)-cx).^2+(y(i)-cy).^2);
%        gx=2*(D-1.5)/D*(cx-x(i));
%        gy=2*(D-1.5)/D*(cy-y(i));
%        GradX=GradX+1/2* gx/(sigma^2 +(D-1.5)^2);
%        GradY=GradY+ 1/2* gy/(sigma^2 +(D-1.5)^2);
%     end
%     G = [GradX, GradY].';
% 
% end