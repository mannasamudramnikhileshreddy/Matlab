
%% References: 

%%          Dolcetti, G., and Krynkin, A. (2020), Matlab Codes for Two-Dimensional
%%      Scattering Surface Reconstruction using Broadband Data, Zenodo

%%          Dolcetti, G., Alkmim, M., Cuenca, J., De Ryck, L., Krynkin, A. (2020), 
%%      Robust Surface Shape Inversion with a Linear Microphones Array, submitted to 
%%      Journal of Sound and Vibrations

%% Author: Giulio Dolcetti g.dolcetti@sheffield.ac.uk
%% Date: 17-07-2020
%% Distribution: Creative Commons Attribution 4.0 International

%%-------------------------------------------------------------------------


%% KirchhoffScattering2D:
%% calculates scattering acoustic potential based on Kirchhoff approximation
%% uses stationary phase approximation along y (surface constant along y)

%% INPUT:

%% CoordMic:        matrix of microphones co-ordinates (Nm x 3)
%%     column1: (x1, x2,..., xM ..., xNm)T   x-co-ordinate  (m)
%%     column2: (y1, y2,..., yM ..., yNm)T   y-co-ordinate  (m)
%%     column3: (z1, z2,..., zM ..., zNm)T   z-co-ordinate  (m)

%% CoordSource:     matrix of source co-ordinates (1 x 5)
%%     column1: xS       x-co-ordinate               (m)
%%     column2: yS       y-co-ordinate               (m)
%%     column3: zS       z-co-ordinate               (m)
%%     column4: psi      angle from horizontal       (rad)
%%     column5: aS       equivalent source radius    (m)


%% CoordSurface:     matrix of surface co-ordinates (Nx x 3)
%%     column1: xx       x-co-ordinate               (m)
%%     column2: yy       y-co-ordinate               (m)
%%     column3: zz       z-co-ordinate               (m)

%% lambda:          acoustic wavelength (m)

%% OUTPUT

%% U:               acoustic potential (-)   (Nm x 1)



function U = KirchhoffScattering2D(CoordMic,CoordSource,CoordSurface,lambda)

    %% non-dimensionalisation
    
k = 2*pi;   % non-dimensional wavenumber
    
xM=repmat(CoordMic(:,1),[1,size(CoordSurface,1)])/lambda;
yM=repmat(CoordMic(:,2),[1,size(CoordSurface,1)])/lambda;
zM=repmat(CoordMic(:,3),[1,size(CoordSurface,1)])/lambda;

xS=repmat(CoordSource(:,1),[size(CoordMic,1),size(CoordSurface,1)])/lambda;
yS=repmat(CoordSource(:,2),[size(CoordMic,1),size(CoordSurface,1)])/lambda;
zS=repmat(CoordSource(:,3),[size(CoordMic,1),size(CoordSurface,1)])/lambda;
psi = CoordSource(:,4);
aS = CoordSource(:,5)/lambda;

xx=repmat(CoordSurface(:,1)',[size(CoordMic,1),1])/lambda;
yy=repmat(CoordSurface(:,2)',[size(CoordMic,1),1])/lambda;
zz=repmat(CoordSurface(:,3)',[size(CoordMic,1),1])/lambda;

%surface slope, central finite difference
etax2=(zz(:,3:end)-zz(:,1:end-2))/(2*mean(diff(xx(1,:))));
etax2=[etax2(:,1),etax2,etax2(:,end)];
etay2 = 0;
    
Dx = mean(diff(xx(1,:)));
    
% components of vectors R1 (source-surface) and R2 (surface receiver)
R1vx=xx-xS; R1vy=yy-yS; R1vz=zz-zS;
R2vx=xM-xx; R2vy=yM-yy; R2vz=zM-zz;

% modulus of R1, R2
R1m=sqrt(R1vx.^2+R1vy.^2+R1vz.^2);
R2m=sqrt(R2vx.^2+R2vy.^2+R2vz.^2);
            
%% directivity
            
% vector of source axis
SourceAxis=[-cos(psi) 0 -sin(psi)];
            
% projections of the distance vectors with respect to source axis
cos1x=R1vx.*SourceAxis(1); cos1y=R1vy.*SourceAxis(2); cos1z=R1vz.*SourceAxis(3);
      
% angle relative to the source axis
theta1=acos(((cos1x+cos1y+cos1z)./R1m).*((cos1x+cos1y+cos1z)./R1m<1)+(1).*((cos1x+cos1y+cos1z)./R1m>=1));

% Directivity (piston with infinite baffle - Morse & Ingard (1986), Theoretical Acoustics, p.381)
Dir=2*besselj(1,k*aS*sin(theta1))./(k*aS*sin(theta1)); Dir(theta1==0)=1;
    

%% Integrand
XX = exp(1i*k*(R1m + R2m));        
    
A2 = 1i/(4*pi) *Dir ./(R1m.*R2m) .*(1 + 1i./(k*R2m))...
                .*(R2vx.*etax2 + R2vy.*etay2 - R2vz)./R2m;
            
A1 = 1i/(4*pi) *Dir ./(R1m.*R2m) .*(1 + 1i./(k*R1m))...
                .*(-R1vx.*etax2 -R1vy.*etay2 + R1vz)./R1m;

A = (A1 + A2); 
            
%% Stationary phase expansion
B = exp(1i*pi/4) .*sqrt(2*pi/k) .*sqrt(R1m.*R2m./(R1m + R2m)) .*A;
            
F = B.*XX;

%% Integration
U = Dx*sum(F,2);
                    


end



