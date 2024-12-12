% Auralius Manurung
% manurunga@yandex.com

clear all;
close all;
clc;

filename = '2poles.gif';

Im = 1; % The current, in Amps
R = 0.5; % in m , this is the outer radius radius of the stator

fig = figure;
hold on;
xlim([-R-R/10 R+R/10]);
ylim([-R-R/10 R+R/10]);
axis equal;

step = 0.04; % spatial grid resolution, try smaller number!
mu_0 = 4*pi*10^-7;

% 2-pole ac motor
% 2-pole motor has 6 wire position, one wire position to antoher makes a
% 60 degress angle

nrotations = 2;
for wt = 0:360*nrotations
    Ia = Im*sind(wt-0);
    Ia_ = -Im*sind(wt-0);
    Ib = Im*sind(wt-120);
    Ib_ = -Im*sind(wt-120);
    Ic = Im*sind(wt-240);
    Ic_ = -Im*sind(wt-240);
    
    xa = [0 -R 0]';
    xa_= [0 R 0]';
    xb = [R*cosd(30) R*sind(30) 0]';
    xb_= [-R*cosd(30) -R*sind(30) 0]';
    xc = [-R*cosd(30) R*sind(30) 0]';
    xc_= [R*cosd(30) -R*sind(30) 0]';
    
    P = [xa xa_ xb xb_ xc xc_];
    I = [Ia Ia_ Ib Ib_ Ic Ic_];
    
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
        ha   = plot(xa(1), xa(2), 'r');
        ha_  = plot(xa_(1), xa_(2), 'r');
        hb   = plot(xb(1), xb(2), 'r');
        hb_  = plot(xb_(1), xb_(2), 'r');
        hc   = plot(xc(1), xc(2), 'r');
        hc_  = plot(xc_(1), xc_(2), 'r');

        I = {Ia, Ia_, Ib, Ib_, Ic, Ic_};
        H = {ha, ha_, hb, hb_, hc, hc_};
        
        imwrite(imind,cm,filename,'gif', 'Loopcount', inf, 'DelayTime', 0.02);
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
        
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime', 0.02);
    end
end



