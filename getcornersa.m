function [xpts,ypts] = getcornersa(amount)
    %figure(get(axishandle, 'Parent'));
    hold on
    xpts=zeros(1,amount);
    ypts=zeros(1,amount);
    for i=1:amount
       [xcur,ycur]=ginput(1);
       xpts(i)=round(xcur);
       ypts(i)=round(ycur);
       if i>1
           plot([xpts(i-1), xpts(i)],[ypts(i-1), ypts(i)],'b');
       end
       plot(xpts(i),ypts(i),'g*');
    end
    %xpts=xcur;
    %ypts=ycur;
end