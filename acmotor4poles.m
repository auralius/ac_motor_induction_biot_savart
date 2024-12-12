% Auralius Manurung
% manurunga@yandex.com

clear all;
close all;
clc;

filename = '4poles.gif';

Im = 1; % The current, in Amps
R = 0.5; % in m , this is the outer radius radius of the stator

fig = figure;
hold on;
xlim([-R-R/10 R+R/10]);
ylim([-R-R/10 R+R/10]);
axis equal;

step = 0.04;
mu_0 = 4*pi*10^-7;

% 4-pole ac motor
% 4-pole motor has 12 wire position, one wire position to antoher makes a
% 30 degress angle

nrotations = 2;
for wt = 0:360*nrotations    
    Ia1 = Im*sind(wt-0);
    Ia1_ = -Im*sind(wt-0);
    Ib1 = Im*sind(wt-120);
    Ib1_ = -Im*sind(wt-120);
    Ic1 = Im*sind(wt-240);
    Ic1_ = -Im*sind(wt-240);
    
    Ia2 = Im*sind(wt-0);
    Ia2_ = -Im*sind(wt-0);
    Ib2 = Im*sind(wt-120);
    Ib2_ = -Im*sind(wt-120);
    Ic2 = Im*sind(wt-240);
    Ic2_ = -Im*sind(wt-240);
    
    xa1 = [0 -R 0]';
    xa1_= [R 0 0]';
    xb1 = [R*cosd(30) -R*sind(30) 0]';
    xb1_= [R*cosd(60) R*sind(60) 0]';
    xc1 = [R*cosd(30) R*sind(30) 0]';
    xc1_= [-R*cosd(60) R*sind(60) 0]';
    
    
    xa2 = [0 R 0]';
    xa2_= [-R 0 0]';
    xb2 = [-R*cosd(30) R*sind(30) 0]';
    xb2_= [-R*cosd(60) -R*sind(60) 0]';
    xc2 = [-R*cosd(30) -R*sind(30) 0]';
    xc2_= [R*cosd(60) -R*sind(60) 0]';
    
    P = [xa1 xa1_ xb1 xb1_ xc1 xc1_ xa2 xa2_ xb2 xb2_ xc2 xc2_];
    I = [Ia1 Ia1_ Ib1 Ib1_ Ic1 Ic1_ Ia2 Ia2_ Ib2 Ib2_ Ic2 Ic2_];
    
    index = 1;
    
    for x = -R:step:R
        for y = -R:step:R
            
            B = [0;0;0];
            
            for k = 1 : length(P)
                if I(k) ~= 0
                    dL = [0;0;I(k)/abs(I(k))];
                else
                    dL = [0;0;0];
                end
                
                r = [x;y;0] - P(:,k);
                r_norm = norm(r);
                
                if r_norm == 0
                    break;
                end
                
                r_hat = r / r_norm;
                
                dB = abs(I(k)) * mu_0 / (4 * pi * r_norm^2) * cross(dL, r_hat);
                B = B + dB;
                
            end
            
            if norm(B) > 1e-10
                B_hat = B/norm(B);
            else
                B_hat=[0;0;0];
            end
            
            X(index) = x;
            Y(index) = y;
            Z(index) = 0;
            U(index) = B_hat(1);
            V(index) = B_hat(2);
            W(index) = B_hat(3);
            
            index = index + 1;
        end
    end
    
    frame = getframe(fig);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    if wt == 0
        h    = quiver3(X,Y,Z, U, V, W);

        ha1  = plot(xa1(1), xa1(2), 'ro');
        ha1_ = plot(xa1_(1), xa1_(2), 'ro');
        hb1  = plot(xb1(1), xb1(2), 'ro');
        hb1_ = plot(xb1_(1), xb1_(2), 'ro');
        hc1  = plot(xc1(1), xc1(2), 'ro');
        hc1_ = plot(xc1_(1), xc1_(2), 'ro');
        
        ha2  = plot(xa2(1), xa2(2), 'ro');
        ha2_ = plot(xa2_(1), xa2_(2), 'ro');
        hb2  = plot(xb2(1), xb2(2), 'ro');
        hb2_ = plot(xb2_(1), xb2_(2), 'ro');
        hc2  = plot(xc2(1), xc2(2), 'ro');
        hc2_ = plot(xc2_(1), xc2_(2), 'ro');

        I = {Ia1, Ia1_, Ib1, Ib1_, Ic1, Ic1_, Ia2, Ia2_, Ib2, Ib2_, Ic2, Ic2_};
        H = {ha1, ha1_, hb1, hb1_, hc1, hc1_, ha2, ha2_, hb2, hb2_, hc2, hc2_};
        
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime', 0.02);
    else
        set(h,'xdata',X,'ydata',Y,'zdata',Z,'udata',U, 'vdata',V,'wdata',W)

        for n = 1 : length(H)
            if I(n) > 0 
                set(H{n}, 'Marker', 'o')
            else
                set(H{n}, 'Marker', 'x')
            end
        end
        
        drawnow
        %pause(0.1);
        
        imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', 0.02);
    end
end

