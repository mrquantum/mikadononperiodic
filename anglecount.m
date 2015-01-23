function theta=anglecount(springlist,XY,angle)

%e1=[tan(angle) 1]/sqrt(1+(tan(angle))^2);
e1=[0 1];
num=numel(XY)/2;
    j=1;

for i=1:numel(springlist(:,1))
     
    
    one=springlist(i,1);
    two=springlist(i,2);
    wlr=springlist(i,3);
    wud=springlist(i,4);
    
    if wud==0 && wlr==0;
    
    y1=XY(one+num);
    y2=XY(two+num);
    
    y1b=XY(one+num);
    y2b=XY(two+num)+wud;
    
    x1=XY(one)+tan(angle)*y1b;
    x2=XY(two)+wlr+tan(angle)*y2b;
    
    v1=[x2-x1 y2-y1]/sqrt((x2-x1)^2+(y2-y1)^2);
    
   
    crosp=-v1(1);
    
    if crosp<0
        theta(j)=-acos(dot(v1,e1));
    end
    if crosp>=0
        theta(j)=acos(dot(v1,e1));
    end
    
    
    if theta(j)>pi/2 && theta(j)<pi
        theta(j)=theta(j)-pi;
    end
    if theta(j)<-pi/2 && theta(j)>-pi
        theta(j)=theta(j)+pi;
    end
    j=j+1;
    
    
    
    
    end
    
end




end

