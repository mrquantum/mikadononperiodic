function plotsprings4(springlist,XY,angle)

e1x=1;
e1y=0;
e2x=tan(angle);
e2y=1;
num=numel(XY)/2;


for gamma=-1:1
    for delta=-1:1
    for(i=1:numel(springlist(:,1)))
        
        one=springlist(i,1);
        two=springlist(i,2);
        wlr=springlist(i,3);
        wud=springlist(i,4);
        
        
        
        y1=XY(one+num)+delta*e2y+gamma*e1y;
        y2=XY(two+num)+wud+delta*e2y+gamma*e1y;
        
        y1b=XY(one+num);
        y2b=XY(two+num)+wud;
        
        x1=XY(one)+gamma*e1x+delta*e2x+tan(angle)*y1b;
        x2=XY(two)+wlr+gamma*e1x+delta*e2x+tan(angle)*y2b;
        
        
        plot([x1,x2],[y1,y2],'k')
        hold on
    end
    end
    
end

end