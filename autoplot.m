close all; clear all; clc;
set(0,'DefaultFigurePosition',[800 100 560 420]);%left bottom width heigth

% path1='/home/mathijs/projects/Mikado/build/mikado1.txt';
% path2='/home/mathijs/projects/Mikado/build/nodes.txt';
% path3='/home/mathijs/projects/Mikado/build/springs.txt';
% path4='/home/mathijs/projects/Mikado/build/conjpoints.txt';
% path5='/home/mathijs/projects/Mikado/build/springboundaries.txt';
% path6='/home/mathijs/projects/Mikado/build/Energy.txt';

path1='/home/mathijs/projects/Mikado/release/mikado1.txt';
path2='/home/mathijs/projects/Mikado/release/nodes.txt';
path3='/home/mathijs/projects/Mikado/release/springs.txt';
path4='/home/mathijs/projects/Mikado/release/conjpoints.txt';
path5='/home/mathijs/projects/Mikado/release/springboundaries.txt';
path6='/home/mathijs/projects/Mikado/release/Energy.txt';
path7='/home/mathijs/projects/Mikado/release/dEda.txt';
path8='/home/mathijs/projects/Mikado/release/rootalpha.txt';


Mikadodata=importdata(path1);
Nodedata=importdata(path2);
Nodes(:,1)=Nodedata(:,2);
Nodes(:,2)=Nodedata(:,3);

springs=importdata(path3);
springs(:,1)=springs(:,1)+1;
springs(:,2)=springs(:,2)+1;
longsp=max(springs(:,5));
shortsp=min(springs(:,5));

 XY=importdata(path4);
 XY=XY';
 anker=[1 .5 .5];

Energy=importdata(path6);

Estrech=Energy(:,1);    
Ebend=Energy(:,2); 
Etot=Energy(:,3); 
lenGrad=Energy(:,4);
conjsteps=1:numel(Ebend);
 
%make plot of the function dEda(a);
cstep=6;
 

% dEda=importdata(path7)';
% roots=importdata(path8);
% alpha=dEda(:,1);
% dEda1=dEda(:,cstep);

% figure
% for i=2:100
   %scatter(roots(i-1),0);
   %hold on
   %plot(alpha,dEda(:,i))
   %hold on
   %plot([-1 1],[0 0])
   %plot([0 0],[-1 1])
   %axis([-.1 .1 -.1 .1])
   %grid on
   %hold off
   %pause(.5)
% end
% figure


%scatter(alpha,dEda1,10,'fill')
% hold on
% plot([-1 1],[0 0],'k');
% plot([0 0],[-.1 .1],'k');
% axis([-.2 .2 -.1 .1])
% xlabel('\alpha','FontSize',18)
%  ylabel('dEda','FontSize',18)
%  title(strcat('Conjstep ',int2str(cstep)),'FontSize',18)

 %grid on;

figure
scatter(conjsteps,Estrech,10,'fill')
hold on
scatter(conjsteps,Ebend,10,'fill');
hold on
scatter(conjsteps,Etot,10,'fill')
hold on;
plot(conjsteps,.1*lenGrad)
legend('Stetching Energy','Bending Energy','Total Energy','Length gradient')


Energy(end,3)
min(Energy(:,3))
axis([0 numel(Energy(:,3))+10 0 0.2])
figure
%set(0,'DefaultFigurePosition',[200 100 560 420]);%left bottom width heigth

%for i=1:length(XY(1,:))
 %   for i=length(XY(1,:)):length(XY(1,:))
   for i=1:1;
       i=8000
    XXYY=XY(:,i);
    X=XY(1:numel(XY(:,i))/2,i);
    Y=XY(numel(XY(:,i))/2+1:end,i);
   %subplot(2,1,1)
   axis([-.01 1.01 -.01 1.01])
    scatter(mod(X,1),mod(Y,1),15,'s','fill','r')
    hold on
    plotsprings(springs,anker,XXYY)
    plotbox(1,1)
    axis equal;
    hold on;
    title(int2str(i),'Fontsize',18)
    end

%     pause(.1);
%     hold off
%     
%     subplot(2,1,2)
%     scatter(i,Ebend(i),'s','k','fill')
%     hold on
%     scatter(i,Estrech(i),'s','r','fill')
%     hold on
%     axis([0 50 0 1])
%     
    
%
%xlabel('x','Fontsize',25)
%ylabel('y','Fontsize',25)
