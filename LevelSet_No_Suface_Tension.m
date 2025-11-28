clear all; clc; close all;
Lx = 1.0; Ly = 1.0;                                                  % domain size 
gx = 0.0; gy = -100.0; rho1 = 2; rho2 = 1; mu = 0.01; gamma = 10;    % parameters
unorth = 0; usouth = 0; veast = 0; vwest = 0;                        % boundary conditions
rad = 0.15; xc = 0.5; yc = 0.7;                                      % initial drop size and location

time = 0.0; plot_freq = 30;

nx = 256; ny = 256; dx = Lx/nx; dy = Ly/ny; dt = 0.0002;
nstep = 2600; maxit = 200; maxError = 0.001; omg = 1.5; Nf = 100;
eps=1.5*dx;
u=zeros(nx+1,ny+2); ut = u   ; uplot = zeros(nx+1,ny+1);
v=zeros(nx+2,ny+1); vt = u   ; vplot = zeros(nx+1,ny+1);
p=zeros(nx+2,ny+2); tmp1 = p ; tmp2  = p;r=p;  chi = p; dirac=p;F=p;
phi1=p;phi2=p;phi3=p;phi4=p;Lap=p;
phi=zeros(nx+2,ny+2);

xf=zeros(1,Nf+2); yf=zeros(1,Nf+2); 
uf=zeros(1,Nf+2); vf=zeros(1,Nf+2);

xh = linspace(0,Lx,nx+1)         ; yh = linspace(0,Ly,ny+1);                % velocity points
x  = linspace(-dx/2,Lx+dx/2,nx+2); y  = linspace(-dy/2,Ly+dy/2,ny+2);       % pressure points

                                               
 fgx=zeros(nx+2,ny+2); fgy=zeros(nx+2,ny+2);
                                          
for l=1:Nf+2; xf(l) = xc - rad*sin(2*pi*(l-1)/Nf); yf(l) = yc + rad*cos(2*pi*(l-1)/Nf); end   % initialize interface.                                  
for i=2:nx+1%Intializing the signed distance function.
    for j=2:ny+1
        for k=1:Nf+2
            d(k)=sqrt((x(i)-xf(k))^2+(y(j)-yf(k))^2);
        end
        phi(i,j)=min(d);
    end
end
phi(1,:) = phi(2,:); phi(nx+2,:) = phi(nx+1,:); phi(:,1) = phi(:,2); phi(:,ny+2) = phi(:,ny+1);
for i = 2:nx+1; for j = 2:ny+1                   
  if((x(i)-xc)^2+(y(j)-yc)^2 < rad^2); phi(i,j) = -phi(i,j); end; 
end; end
r = zeros(nx+2,ny+2) + rho2;
 for i=1:nx+2;for j=1:ny+2;
         if phi(i,j)<-eps
             chi(i,j)=1;
             r(i,j)=rho1;
         end
         if phi(i,j)>=-eps&&phi(i,j)<eps
             chi(i,j)=0.5+0.5*(phi(i,j)/eps)+(0.5/(pi))*(sin((pi*phi(i,j))/eps));
         end
          if phi(i,j)>eps
             chi(i,j)=0;
            
          end
     end
 end

        
     
for is=1:nstep
    phi_old=phi;
 for i=2:nx+1%Upwind
     for j=2:ny+1
         phi(i,j)=phi(i,j)-dt*(max(u(i,j),0)*(((phi(i,j)-phi(i-1,j)))/(dx))+(min(u(i,j),0)*(phi(i+1,j)-phi(i,j))/(dx))...
                  +(max(v(i,j),0)*(phi(i,j)-phi(i,j-1))/(dy))+(min(v(i,j),0)*(phi(i,j+1)-phi(i,j))/(dy)));
     end
 end

phi(1,:) = phi(2,:); phi(nx+2,:) = phi(nx+1,:); phi(:,1) = phi(:,2); phi(:,ny+2) = phi(:,ny+1); 
 
 phi0=phi;
 F=phi0./(sqrt(phi0.^2+eps^2));
  

 for tim=1:10
     for i=2:nx+1
         for j=2:ny+1%Solving Hamilton Jacobi,Reinitalization
             phi(i,j)=phi(i,j)-dt*(((max(F(i,j),0))*sqrt(((max((phi(i,j)-phi(i-1,j))/(dx),0)^2))+((min((phi(i+1,j)-phi(i,j))/(dx),0)^2)) ...
                 +((max((phi(i,j)-phi(i,j-1))/(dy),0))^2)+((min((phi(i,j+1)-phi(i,j))/(dy),0))^2)))...
                 +((min(F(i,j),0))*(sqrt((max((phi(i+1,j)-phi(i,j))/(dx),0)^2)+((min((phi(i,j)-phi(i-1,j))/(dx),0))^2)...
                 +((max((phi(i,j+1)-phi(i,j))/(dy),0))^2)+((min((phi(i,j)-phi(i,j-1))/(dy),0)^2))))))+dt*(F(i,j));
                 
         end
     end
 phi(1,:) = phi(2,:); phi(nx+2,:) = phi(nx+1,:); phi(:,1) = phi(:,2); phi(:,ny+2) = phi(:,ny+1);
    
 end
 for i=2:nx+1
     for j=2:ny+1
         delphi(i,j)=sqrt(((phi(i+1,j)-phi(i,j))/(dx))^2+((phi(i,j+1)-phi(i,j))/(dy))^2);
     end
 end

 
   for i=1:nx+2;for j=1:ny+2; %Smeared out characteristic function
         if phi(i,j)<-eps
             chi(i,j)=1;
         end
         if phi(i,j)>=-eps&phi(i,j)<=eps
             chi(i,j)=0.5+0.5*(phi(i,j)/eps)+(0.5/(pi))*(sin((pi*phi(i,j))/eps));
         end
          if phi(i,j)>eps
             chi(i,j)=0;
          end
     end
   end
  ro = r;
  r  = rho1*chi + rho2*(1-chi);  % obtain density from charact func
 for i=1:nx+2;
      for j=1:ny+2;  %Dirac Delta
          if phi(i,j)<-eps
              dirac(i,j)=0;
          end
          if phi(i,j)>=-eps&&phi(i,j)<=eps
              dirac(i,j)=(0.5/eps)+(0.5/eps)*(cos((pi*phi(i,j))/eps));
          end
          if phi(i,j)>eps
              dirac(i,j)=0;
          end
      end
 end
 for i=2:nx+1%Laplacian of phi
     for j=2:ny+1
         Lap(i,j)=(((phi(i+1,j)+phi(i-1,j)-2*phi(i,j))/(dx^2))+((phi(i,j+1)+phi(i,j-1)-2*phi(i,j))/(dy^2)));%Laplacian
     end
 end

%  for i=2:nx+1;  %Surface Tension
%      for j=2:ny+1;
%          fgx(i,j)=gamma*dirac(i,j)*Lap(i,j)*((phi(i+1,j)-phi(i-1,j))/(2*dx));
%                   
%          fgy(i,j)=gamma*dirac(i,j)*Lap(i,j)*((phi(i,j+1)-phi(i,j-1))/(2*dy));    
%              
%                   
%      end
%  end

  fgx(1:nx+2,2) = fgx(1:nx+2,2) + fgx(1:nx+2,1); fgx(1:nx+2,ny+1) = fgx(1:nx+2,ny+1) + fgx(1:nx+2,ny+2);  % bring all forces to interior
  fgy(2,1:ny+2) = fgy(2,1:ny+2) + fgy(1,1:ny+2); fgy(nx+1,1:ny+2) = fgy(nx+1,1:ny+2) + fgy(nx+2,1:ny+2);  
 
   u(1:nx+1,1) = 2*usouth-u(1:nx+1,2); u(1:nx+1,ny+2) = 2*unorth-u(1:nx+1,ny+1); % tangential vel BC
   v(1,1:ny+1) = 2*vwest -v(2,1:ny+1); v(nx+2,1:ny+1) = 2*veast -v(nx+1,1:ny+1); % tangential vel BC
 

  for i=2:nx; for j=2:ny+1     % temporary u-velocity (boundary values are not touched)%  %Adding surface tension to ut and vt
    ut(i,j) = (2.0/(r(i+1,j)+r(i,j)))*(0.5*(ro(i+1,j)+ro(i,j))*u(i,j)+ dt* (...
            - (0.25/dx)*(ro(i+1,j)*(u(i+1,j)+u(i,j))^2-ro(i,j)*(u(i,j)+u(i-1,j))^2)...
            - (0.0625/dy)*( (ro(i,j)+ro(i+1,j)+ro(i,j+1)+ro(i+1,j+1))*(u(i,j+1)+u(i,j))*(v(i+1,j)+v(i,j)) ...
            - (ro(i,j)+ro(i+1,j)+ro(i+1,j-1)+ro(i,j-1))*(u(i,j)+u(i,j-1))*(v(i+1,j-1)+v(i,j-1)))...
            + mu*((u(i+1,j)-2*u(i,j)+u(i-1,j))/dx^2+ (u(i,j+1)-2*u(i,j)+u(i,j-1))/dy^2)...
            + 0.5*(ro(i+1,j)+ro(i,j))*gx +fgx(i,j)   ) ); 
  end; end

  for i=2:nx+1; for j=2:ny       % temporary v-velocity (boundary values are not touched)
    vt(i,j) = (2.0/(r(i,j+1)+r(i,j)))*(0.5*(ro(i,j+1)+ro(i,j))*v(i,j)+ dt* (...     
            - (0.0625/dx)*( (ro(i,j)+ro(i+1,j)+ro(i+1,j+1)+ro(i,j+1))*(u(i,j)+u(i,j+1))*(v(i,j)+v(i+1,j)) ...
            - (ro(i,j)+ro(i,j+1)+ro(i-1,j+1)+ro(i-1,j))*(u(i-1,j+1)+u(i-1,j))*(v(i,j)+v(i-1,j)) )...                                 
            - (0.25/dy)*(ro(i,j+1)*(v(i,j+1)+v(i,j))^2-ro(i,j)*(v(i,j)+v(i,j-1))^2 )...
            + mu*((v(i+1,j)-2*v(i,j)+v(i-1,j))/dx^2+(v(i,j+1)-2*v(i,j)+v(i,j-1))/dy^2)...
            + 0.5*(ro(i,j+1)+ro(i,j))*gy +fgy(i,j) ) );    
  end; end     

  for i = 2:nx+1; for j = 2:ny+1
    tmp1(i,j) = (0.5/dt)*( (ut(i,j)-ut(i-1,j))/dx+(vt(i,j)-vt(i,j-1))/dy );
    tmp2(i,j) =1/( (1/dx)*(1/(dx*(r(i+1,j)+r(i,j)))+ 1/(dx*(r(i-1,j)+r(i,j))) )+ ...
                   (1/dy)*(1/(dy*(r(i,j+1)+r(i,j)))+ 1/(dy*(r(i,j-1)+r(i,j))) )   );
  end; end

  for it = 1:maxit	               % solve for pressure by SOR
    pold   = p;
    p(1,:) = p(2,:); p(nx+2,:) = p(nx+1,:); p(:,1) = p(:,2); p(:,ny+2) = p(:,ny+1); % set gosht values
    for i=2:nx+1; for j=2:ny+1
      p(i,j) = (1.0-omg)*p(i,j) + omg*tmp2(i,j)*(        ...
      (1/dx)*( p(i+1,j)/(dx*(r(i+1,j)+r(i,j)))+ p(i-1,j)/(dx*(r(i-1,j)+r(i,j))) )+    ...
      (1/dy)*( p(i,j+1)/(dy*(r(i,j+1)+r(i,j)))+ p(i,j-1)/(dy*(r(i,j-1)+r(i,j))) ) - tmp1(i,j));
    end; end
    if max(max(abs(pold-p))) < maxError; break;  end
  end
                                      
  for i=2:nx; for j=2:ny+1   % correct the u-velocity 
    u(i,j)=ut(i,j)-dt*(2.0/dx)*(p(i+1,j)-p(i,j))/(r(i+1,j)+r(i,j));
  end; end
      
  for i=2:nx+1; for j=2:ny   % correct the v-velocity 
    v(i,j)=vt(i,j)-dt*(2.0/dy)*(p(i,j+1)-p(i,j))/(r(i,j+1)+r(i,j));
  end; end
[r1,c1]=find(phi>-eps&phi<eps);
 xi=x(r1);
 yi=y(c1);

  time = time+dt;                       
  if (mod(is,plot_freq)==0) | (is==1);                         % plot solution
  uu(1:nx+1,1:ny+1)=0.5*(u(1:nx+1,2:ny+2)+u(1:nx+1,1:ny+1));
  vv(1:nx+1,1:ny+1)=0.5*(v(2:nx+2,1:ny+1)+v(1:nx+1,1:ny+1));
 figure(1); contour(x,y,r'); axis equal;title(sprintf('Geometry t=%0.2f',time)); axis([0 Lx 0 Ly]); hold on;quiver(xh,yh,uu',vv','r'); hold on 
 plot(xi,yi,'k.','linewidth',3);xlim([0,1]);ylim([0,1]); drawnow; hold off;
 
  %figure(2); plot(xf(1:Nf)); hold off;
  end

end                  
