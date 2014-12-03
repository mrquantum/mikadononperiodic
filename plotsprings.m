function []=plotsprings(spring,anker,coords)
s=spring;
%s(:,5)=ones(numel(s(1,:),1));
x=coords(1:length(coords)/2); y=coords(1+length(coords)/2:end);

longsp=0.01;
shortsp=0.0001;
%s(:,5)

for i=1:length(s(:,1))
        if s(i,3)==0 && s(i,4)==0
            %plot([mod(x(s(i,1)),1),mod(x(s(i,2)),1)],[y(s(i,1)),y(s(i,2))],'k','LineWidth',1)
            if s(i,5)==longsp
            plot([x(s(i,1)),x(s(i,2))],[y(s(i,1)),y(s(i,2))],'b','LineWidth',1)
            end           
            if s(i,5)==shortsp
            plot([x(s(i,1)),x(s(i,2))],[y(s(i,1)),y(s(i,2))],'k','LineWidth',2)
            end
            
            hold on
        end
        
        if abs(s(i,3))==1 && s(i,4)==0
           if s(:,5)==longsp
            plot([x(s(i,1)),x(s(i,2))+s(i,3)*1],[y(s(i,1)),y(s(i,2))],'b','LineWidth',1)
            hold on
            plot([x(s(i,1))-s(i,3)*1,x(s(i,2))],[y(s(i,1)),y(s(i,2))],'b','LineWidth',1)
           end
            if s(:,5)==shortsp
            plot([x(s(i,1)),x(s(i,2))+s(i,3)*1],[y(s(i,1)),y(s(i,2))],'k','LineWidth',2)
            hold on
            plot([x(s(i,1))-s(i,3)*1,x(s(i,2))],[y(s(i,1)),y(s(i,2))],'k','LineWidth',2)
           end
        end
        
        if s(i,3)==0 && abs(s(i,4))==1
            if s(i,5)==longsp
            plot([x(s(i,1)),x(s(i,2))],[y(s(i,1)),y(s(i,2))+s(i,4)*1],'b','LineWidth',1)
            hold on
            plot([x(s(i,1)),x(s(i,2))],[y(s(i,1))-s(i,4)*1,y(s(i,2))],'b','LineWidth',1)
            end
            if s(i,5)==shortsp
                 plot([x(s(i,1)),x(s(i,2))],[y(s(i,1)),y(s(i,2))+s(i,4)*1],'k','LineWidth',2)
            hold on
            plot([x(s(i,1)),x(s(i,2))],[y(s(i,1))-s(i,4)*1,y(s(i,2))],'k','LineWidth',2)
            end
                
            end
        
        
end