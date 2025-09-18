function theta_1x = Protocol_Wm(deltax,a,V0,offset)
% This fucnction computes the protocol with a measure uncertanity \Delta x
L=2;
xdis=linspace(0,L,1000); % x is shifted to the interval (0,L=2)
theta_1x=zeros(1, length(xdis));
p=0;
%Determinamos xl y xr
xl=a*(V0-offset)/V0;
xr=(offset*(L-a)+V0*a)/V0;
for i=1:length(xdis)
    x=xdis(i);
    l_u= mod(x+deltax,L); % upper limit of the measure interval
    l_l=mod(x-deltax,L);  % lower limit of the measure interval
    if (deltax==0)
        p=(xl<x)&&(x<xr);
        theta_1x(i)=p;
    else
        if (l_l<=l_u)
            if (l_u<=xl)
                p=0;
            elseif (xl<=l_u)&&(l_u<=xr)
                if (l_l<=xl)
                    p=(l_u-xl)/(2*deltax);
                else %l_l>xl
                    p=(l_u-l_l)/(2*deltax);
                end
            else %l_u>xr
                if (l_l<=xl)
                    p=(xr-xl)/(2*deltax);
                elseif (xl<=l_l)&&(l_l<=xr)
                    p=(xr-l_l)/(2*deltax);
                else %l_l>xr
                    p=0;
                end
            end
        else %l_l>l_u
            if (l_l<=xl)
                p=(xr-xl)/(2*deltax);
            elseif (xl<=l_l)&&(l_l<=xr)
                if (l_u<=xl)
                    p=(xr-l_l)/(2*deltax);
                else %l_u>xl
                    p= (xr-xl+l_u-l_l)/(2*deltax);
                end
            else %l_l>xr
                if (l_u<=xl)
                    p=0;
                elseif (xl<=l_u)&&(l_u<=xr)
                    p=(l_u-xl)/(2*deltax);
                else %l_u>xr
                    p=(xr-xl)/(2*deltax);
                end
            end
        end
        theta_1x(i)=p;
    end
end
end


    