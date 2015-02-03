close all; clear all; clc;
set(0,'DefaultFigurePosition',[800 100 560 420]);%left bottom width heigth

%Path to different text-files
path1='/home/mathijs/projects/MikadoSlanted/release/data/mikado.txt';
path2='/home/mathijs/projects/MikadoSlanted/release/data/nodes.txt';
path3='/home/mathijs/projects/MikadoSlanted/release/data/springs.txt';
path4='/home/mathijs/projects/MikadoSlanted/release/data/conjpoints.txt';
path5='/home/mathijs/projects/MikadoSlanted/release/data/springboundaries.txt';
path6='/home/mathijs/projects/MikadoSlanted/release/data/Energy.txt';
path7='/home/mathijs/projects/MikadoSlanted/release/data/dEda.txt';
path8='/home/mathijs/projects/MikadoSlanted/release/data/rootalpha.txt';
path9='/home/mathijs/projects/MikadoSlanted/release/data/shearcoordinates.txt';
path10='/home/mathijs/projects/MikadoSlanted/release/data/shearenergy.txt';

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

Energy=importdata(path6);
shearedpos=importdata(path9)';
shearenergy=importdata(path10);


%%Energy as function of shearingangle
figure
plot(tan(shearenergy(:,2)),shearenergy(:,1));
xlabel('\gamma=tan(\phi)','FontSize',18)
ylabel('Energy','FontSize',18)

Estrech=Energy(:,1);
Ebend=Energy(:,2);
Etot=Energy(:,3);
lenGrad=Energy(:,4);
conjsteps=1:numel(Ebend);
grid on
axis([-.4 .4 .115 .118])


%Here the energy is plotted as a function of the conjugate steps
% figure
% semilogy(conjsteps(1:150),Estrech(151:300),'r','LineWidth',2)
% hold on
% semilogy(conjsteps(1:150),Ebend(151:300),'k','LineWidth',2);
% hold on
% semilogy(conjsteps(1:150),Etot(151:300),'b','LineWidth',2)
% hold on;
% %semilogy(conjsteps,.1*lenGrad,'LineWidth',2)
% legend('Stetching Energy','Bending Energy','Total Energy')%,'Length gradient')
% %axis([0 numel(Energy(:,3))+10 1e-3 .1])
% grid on
% xlabel('Conjugate steps','FontSize',18)
% ylabel('Energy','FontSize',18)
% 
% 
% %Here the network is plotted (Nodes and springs)
% %phis=0:.05:0.45;
% phis=0:0.01:0.24;
% phis2=.24:-.01:-.24;
% phis=horzcat(phis,phis2);
% figure
% i=1;
% XXYY=shearedpos(:,i);
% 
% %%Now We plot everything in normal coordinates.
% phi=phis(i);
% 
% e1y=0;
% e2x=tan(phi);
% e2y=1;
% e1x=1;
% plotsprings4(springs,XXYY,phi)
% hold on
% plot([0 e1x],[0 e1y],'LineWidth',2); plot([0,e2x],[0,e2y],'LineWidth',2)
% plot([1 1+e2x],[0,e2y],'LineWidth',2);plot([tan(phi) 1+tan(phi)],[1 1],'LineWidth',2)
% plot([-1 2],[0,0]);plot([-1,2],[-1,-1]);
% plot([-1,2],[1,1])
% plot([0,0],[-1,2])
% plot([1,1],[-1,2])
% axis([-1 2 -1 2])
% grid on
% 
% 
% 
% axis equal;
% hold on;
% title(int2str(i),'Fontsize',18)
% xlabel('x','Fontsize',25)
% ylabel('y','Fontsize',25)
% title(strcat('\phi=',num2str(phi)),'FontSize',18)
% 
% 
% 
% 
% 
% % Frame(i)=getframe(gcf)
% % clf
% % end
% % movie2avi(Frame,'shearing.avi','fps',1)
% 
% 
% % th=anglecount(springs,XXYY,phi);
% % figure
% % hist(th,25)
% % title(strcat('For tilted under \phi=',num2str(phi)),'FontSize',18)
% % xlabel('\theta','FontSize',18)
% % ylabel('Count','FontSize',18)
% 
% 
