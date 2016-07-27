%fdtd analysis
%
%
%fdtd engine
close all;
clc;
clear all;

%units
 meters = 1;
centimeters = 1e-2 *  meters;
millimeters = 1e-3 * meters;
inches = 2.54 * centimeters;
feet = 12 * inches;

seconds  = 1;
hertz     = 1/seconds;
kilohertz = 1e-3 * hertz;
megahertz = 1e-6 * hertz;
gigahertz = 1e-9 * hertz;

%constant

c0 = 299792458 * meters/seconds;
e0 = 8.8541878176e-12 * 1/meters;
u0 = 1.2566370614e-6 * 1/meters;

%open a figure window
figure('Color','w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%dashboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Source parameters
fmax = 5.0 * gigahertz;

%grid parameters

nmax = 1;
NLAM = 10;
NBUFZ = [100,100];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMPUTE OPTIMISED GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%NOMINAL RESOLUTION

lam0 = c0/fmax;
dz = lam0/nmax/NLAM;

%compute grid size

Nz = sum(NBUFZ) + 3;


%COMPUTE GRID AXIS
za = [0: Nz-1]*dz;

%dz = 0.006 * meters;
%Nz = 200;
%dt =1e-11*seconds;
%STEPS = 1000;


% BUILD DEVICE ON THE GRID


% INITIALIZE MATERIALS TO FREE SPACE

 ER = ones(1,Nz);
 UR = ones(1,Nz);
 
 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%compute the source
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 %compute the time step(dt)

 nbc = sqrt(ER(1)*UR(1));
 dt  = nbc*dz/(2*c0);
 
 
 %compute the source parameters
 
 tau = 0.5/fmax;
 t0 = 5* tau;
 
 
 %compute number of time steps
 
 tprop = nmax*Nz*dz/c0;
 t = 2*t0 + 3*tprop;
 STEPS = ceil(t/dt);
 
 
 % compute the source
 
 t = [0:STEPS-1]*dt;
nz_src = round(Nz/2);
Esrc = exp(-((t-t0)/tau).^2);

%plot(t,Esrc);
%break;
 
 
 
 
 %COMPUTE UPDATE COEFFICIENTS
 
 mEy = (c0*dt)./ER;
 mHx = (c0*dt)./UR;
 
 %initialize fields
 
 Ey = zeros(1,Nz);
 Hx = zeros(1,Nz);
 
 
 %initialize boundry conditions
 H1=0; H2=0; H3=0;
 E1=0; E2=0; E3=0;
 
 
 
 
 
 % PERFORM FDTD ANALYSIS
 
 for T = 1:STEPS
     
     %UPDTATE H FROM E
     
     for nz = 1 : Nz-1
         
         Hx(nz) =  Hx(nz) + mHx(nz)*( Ey(nz+1)-Ey(nz) )/dz;
     end
      Hx(Nz) =  Hx(Nz) + mHx(Nz)*( E3 -Ey(Nz) )/dz
      
      
      %REcord H-field at the boundry
      
      H3=H2;
      H2=H1;
      H1=Hx(1);
      
      
      % UPDATE E FROM H
      
      Ey(1) = Ey(1) +  mEy(1)* ( Hx(1) - H3 )/dz;
      for nz = 2: Nz
          Ey(nz) = Ey(nz) +  mEy(nz)*( Hx(nz) - Hx(nz-1) )/dz;
      end
      
      %RECORD EFFIELD AT THE BOUNDRY
      
      E3=E2;
      E2=E1;
      E1=Ey(Nz);
      
      
 end
 
 
 %inject Esource
 
 Ey(nz_src) =   Ey(nz_src) +  Esrc(T); 
 
 
 
      % SHOW status
      if ~mod(T,1)
          
          %show fields
          draw1d(ER,Ey,Hx,dz);
          xlim([dz,Nz*dz]);
          xlabel('z');
          title(['field at step' num2str(T) 'of' num2str(STEPS)']);
          
          
          % draw graphics
          drawnow;
      end
      
         
         
         
     