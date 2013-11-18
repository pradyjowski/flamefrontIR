function [xpts,ypts] = getpointsa()

    %figure(get(axishandle, 'Parent'));
    hold on
    xpts=[];    ypts=[];    splineh=[];
    i=1;
    [xcur,ycur,button]=ginput(1);
    plot(xcur,ycur,'g*');
    while button~=3
       xpts(i)=round(xcur);
       ypts(i)=round(ycur);
       i=i+1;
       [xcur,ycur,button]=ginput(1);
       if button==3 
           plot(xpts(i-1),ypts(i-1),'g*');
       else
           plot(xcur,ycur,'go');
           %plot([xcur, xpts(i-1)],[ycur, ypts(i-1)],'b');
           t = 1:i;
           ts = 1: 0.25 : i;
           xs = spline(t, [xpts, xcur], ts);
           ys = spline(t, [ypts, ycur], ts);
           if ~ isempty(splineh)
               delete(splineh);
           end;  
           splineh=plot(xs,ys,'b');
           
       end
       
    end
    
    xpts(i:end)=[];
    ypts(i:end)=[];
    %xpts=xcur;
    %ypts=ycur;
end