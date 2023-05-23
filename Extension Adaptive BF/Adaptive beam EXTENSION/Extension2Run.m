clc
clear all
close all

%   Definitions:
ensemble    = 100;                          % number of realizations within the ensemble
K           = 1000;                          % number of iterations
H           = [0.32+0.21*1i,-0.3+0.7*1i,0.5-0.8*1i,0.2+0.5*1i].';
Wo          = H;                            % unknown system
sigma_n2    = 0.04;                         % noise power
N           = 4;                            % number of coefficients of the adaptive filter
mu         =[0.1];                           % convergence factor (step)  (0 < mu < 1)
gamma       = 1e-12;                        % small positive constant to avoid singularity
%   Initializing & Allocating memory:
W       = ones(N,(K+1),ensemble);   % coefficient vector for each iteration and realization; 

%   Computing:
for l=1:ensemble

    X       = zeros(N,1);               % input at a certain iteration (tapped delay line)
    d       = zeros(1,K);               % desired signal

    x        = (sign(randn(K,1)) + j*sign(randn(K,1)))./sqrt(2); % Creating the input signal (normalized)
    sigma_x2 = var(x);                                           % signal power = 1
    n        = sqrt(sigma_n2/2)*(randn(K,1)+j*randn(K,1));       % complex noise

    for k=1:K

        X       =   [x(k,1)
                     X(1:(N-1),1)];              % input signal (tapped delay line)

        d(k)    =   (Wo'*X(:,1))+n(k);           % desired signal

    end

    S   =   struct('step',mu,'filterOrderNo',(N-1),'initialCoefficients',...
                   W(:,1,l),'gamma',gamma);
    [y,e,W(:,:,l)]  =   MODIFIEDLMS(d,transpose(x),S);

end
%   Averaging:
W_av = sum(W,3)/ensemble;
%   Plotting:
figure,
plot(real(W_av(1,:)),'DisplayName','mu=0.1');
title('Mold of vector w');
xlabel('iterations, K'); ylabel('Mold');
hold on
mu         =[0.05];                           % convergence factor (step)  (0 < mu < 1)

%   Initializing & Allocating memory:
W1       = ones(N,(K+1),ensemble);   % coefficient vector for each iteration and realization;


%   Computing:
for l=1:ensemble

    X       = zeros(N,1);               % input at a certain iteration (tapped delay line)
    d       = zeros(1,K);               % desired signal

    x        = (sign(randn(K,1)) + 1i*sign(randn(K,1)))./sqrt(2); % Creating the input signal (normalized)
    sigma_x2 = var(x);                                           % signal power = 1
    n        = sqrt(sigma_n2/2)*(randn(K,1)+1i*randn(K,1));       % complex noise

    for k=1:K

        X       =   [x(k,1)
                     X(1:(N-1),1)];              % input signal (tapped delay line)

        d(k)    =   (Wo'*X(:,1))+n(k);           % desired signal

    end

    S   =   struct('step',mu,'filterOrderNo',(N-1),'initialCoefficients',...
                   W1(:,1,l),'gamma',gamma);
    [y,e1,W1(:,:,l)]  =  MODIFIEDLMS(d,transpose(x),S);
end
%   Averaging:
W_av1 = sum(W1,3)/ensemble;
plot(real(W_av1(1,:)),'DisplayName','mu=0.05'),...
title('Mold of vector w');
xlabel('iterations, K'); ylabel('Mold');

hold on
mu         =[0.02];                           % convergence factor (step)  (0 < mu < 1)

%   Initializing & Allocating memory:
W2       = ones(N,(K+1),ensemble);   % coefficient vector for each iteration and realization; 

%   Computing:
for l=1:ensemble

    X       = zeros(N,1);               % input at a certain iteration (tapped delay line)
    d       = zeros(1,K);               % desired signal

    x        = (sign(randn(K,1)) + 1i*sign(randn(K,1)))./sqrt(2); % Creating the input signal (normalized)
    sigma_x2 = var(x);                                           % signal power = 1
    n        = sqrt(sigma_n2/2)*(randn(K,1)+1i*randn(K,1));       % complex noise

    for k=1:K

        X       =   [x(k,1)
                     X(1:(N-1),1)];              % input signal (tapped delay line)

        d(k)    =   (Wo'*X(:,1))+n(k);           % desired signal

    end

    S   =   struct('step',mu,'filterOrderNo',(N-1),'initialCoefficients',...
                   W2(:,1,l),'gamma',gamma);
    [y2,e2,W2(:,:,l)]  =   MODIFIEDLMS(d,transpose(x),S);

end
%   Averaging:
W_av1 = sum(W2,3)/ensemble;
hold on
plot(real(W_av1(1,:)),'DisplayName','mu=0.02'),...
legend
title('Mold of vector w');
xlabel('iterations, K'); ylabel('Mold');

                          

W22       = ones(N,(K+1),ensemble);   % coefficient vector for each iteration and realization;
wo         =[0.3];
%   Computing:
for l=1:ensemble

    X       = zeros(N,1);               % input at a certain iteration (tapped delay line)
    d       = zeros(1,K);               % desired signal

    x        = (sign(randn(K,1)) + 1i*sign(randn(K,1)))./sqrt(2); % Creating the input signal (normalized)
    sigma_x2 = var(x);                                           % signal power = 1
    n        = sqrt(sigma_n2/2)*(randn(K,1)+1i*randn(K,1));       % complex noise

    for k=1:K

        X       =   [x(k,1)
                     X(1:(N-1),1)];              % input signal (tapped delay line)

        d(k)    =   (Wo'*X(:,1))+n(k);           % desired signal

    end

    S   =   struct('step',wo,'filterOrderNo',(N-1),'initialCoefficients',...
                   W22(:,1,l),'gamma',gamma);
    [y2,e2,W22(:,:,l)]  =   MODIFIEDLMS(d,transpose(x),S);

end
%   Averaging:
W_av11 = sum(W22,3)/ensemble;
hold on
plot(real(W_av11(1,:)),'DisplayName','Optimal Solution'),...
legend
title('Mold of vector w');
xlabel('iterations, K'); ylabel('Mold');
hold off
%channel system order
M = 4 ;
N=1000;
inp = randn(N,1);
n = randn(N,1);
[b,a] = butter(2,0.25);                             %Create a low pass filter of second order with cutoff frequency 0.25
Gz = tf(b,a,.1);                                    %Create discrete-time transfer function with undetermined sample time                                                   %and numerator b and a, which are the transfer function coefficients of the                                                    %2nd order butterworth filter with cutoff freq 0.25
                                                    % if you use ldiv this will give h :filter weights to be
                                        %This is the actual filter that we are trying to recreate
h=[0.32+0.21*1i,-0.3+0.7*1i,0.5-0.8*1i,0.2+0.5*1i];
y = lsim(Gz,inp);                                   %This simulates the time response of the system Gz given random input inp
                                                    %add some noise
n = n * std(y)/(10*std(n));
d = y + n;
totallength=size(d,1);
                                                    %Take 60 points for training
N=60 ;
                                                    %begin of algorithm
w = zeros ( M  , 1 ) ;
for n = M : N 
	u = inp(n:-1:n-M+1) ;
    y(n)= w' * u;                                   %Compute y(n) with the weights and the input
    e(n) = d(n) - y(n) ;                            %Compute the error: d(n) is the true input
                                                    %Start with big mu for speeding the convergence then slow down
                                                    %to reach the correct weights
    if n < 20
        mu=0.05;
    else
        mu=0.15;
    end
                                                    %The update corresponds to the LMS update rule, which can be derived
                                                    %to reach the correct weights
                                                    %using the gradient of the cost function and the approximation of the
                                                                                           %over n points. 
	w = w + mu * u * e(n) ;                         
end 
                                                    %Check of results on
                                                    %the rest of the data
for n =  N+1 : 63
    u = inp(n:-1:n-M+1);
    y(n) = w' * u ; 
	msg = round(y(n));
    modmsg = qammod(msg,4,'bin');
    input = 0:1;
    const = qammod(input,4,'bin');
    tblen = 50; chcffs = [100 ;100.5];
    filtmsg = filter(chcffs,1,modmsg);   
    [rx1, ~, ~, ~] = mlseeq(filtmsg,chcffs,const,tblen,'cont');
    y(n) = isequal(rx1(tblen+1:end),modmsg(1:end-tblen));

                                    %Compute the results with the obtained weight vector
    
    e(n) = d(n) - y(n) ;                            %Compute the error
end
figure
semilogy(sort(abs(e),'descend')) ;                                %Create semi log scale plot (only y axis)
load e1
hold on
semilogy(sort(abs(e1),'descend')) ;                                %Create semi log scale plot (only y axis)
title('Mean Square Error (MSE)') ;
xlabel('Iterations')
ylabel('MSE')
legend('Proposed','Existing')
NumofAntenna = 3; % Number of antennas in the array
NumofSamples = 100; % Number of bits to be transmitted
SigmaSystem = 0.1; % System Noise Variance
theta_x = 2* (pi/180); % direction of signal x
theta_n1 = 40 * (pi/180); % direction of noise source 1
theta_n2 = -40 * (pi/180); % direction of noise source 2
% TIME SETTINGS
theta = pi*[-1:0.005:1];
BitRate = 100;
SimFreq = 4*BitRate; % Simulation frequency
Ts = 1/SimFreq; % Simulation sample period
for k=1:NumofSamples
q=randperm(2);
Data(k)=-1^q(1);
end
Data = upsample(Data, SimFreq/BitRate); % Upsample data
t = Ts:Ts:(length(Data)/SimFreq); % Timeline
faz=(cumsum(Data))/8;
signal_x = cos(pi*faz)+j*sin(pi*faz); % The signal to be received
% GENERATE INTERFERER NOISE -> uniform phase (-pi,pi), gaussian amplitude
% distribution(magnitude 1)
signal_n1 = normrnd(0,1,1,length(t)).*exp (j*(unifrnd(-pi,pi,1,length(t))));
signal_n2 = normrnd(0,1,1,length(t)).*exp (j*(unifrnd(-pi,pi,1,length(t))));
% GENERATE SYSTEM NOISES for EACH ANTENNA -> uniform phase (-pi,pi),
% gaussian
% amplitude distribution(magnitude 1)
noise = zeros(NumofAntenna, length(t));
for i = 0:NumofAntenna-1,
noise(i+1,:) = normrnd(0,SigmaSystem,1,length(t)).*exp (j*(unifrnd(-pi,pi,1,length(t))));
end;
% ARRAY RESPONSES for DESIRED SIGNAL (X) and INTERFERER NOISES (N1
% and N2)
Kd = pi; % It is assumed that antennas are seperated by lambda/2.
response_x = zeros(1,NumofAntenna);
response_n1 = zeros(1,NumofAntenna);
response_n2 = zeros(1,NumofAntenna);
for k = 0:NumofAntenna-1
    response_x(k+1) = exp(j*k*Kd*sin(theta_x));
response_n1(k+1) = exp(j*k*Kd*sin(theta_n1));
response_n2(k+1) = exp(j*k*Kd*sin(theta_n2));
end;
% TOTAL RECEIVED SIGNAL (SUM of X.*Hx, N1.*Hn1 and N2.*Hn2)
x24 = zeros(NumofAntenna, length(t));
n1 = zeros(NumofAntenna, length(t));
n2 = zeros(NumofAntenna, length(t));
for i = 0:NumofAntenna-1,
x24(i+1,:) = signal_x .* response_x(i+1); % received signal from signal source x
n1(i+1,:) = signal_n1 .* response_n1(i+1); % received signal from noise source n1
n2(i+1,:) = signal_n2 .* response_n2(i+1); % received signal from noise source n2
end;
signal_ns = (noise + n1+n2+x24); % total received signal

% EVALUATUING WEIGHTs THOSE SATISFY BEAMFORMING at DESIRED
% DIRECTION
y24 = zeros(1,length(t)); % output
mu24 = 0.05; % gradient constant
e24 = zeros(1,length(t)); % error
w24 = zeros(1,NumofAntenna); w24(1)=eps; % weights
for i=0:length(t)-1,
y24(i+1) = w24 * signal_ns(:,i+1);
e24(i+1) = y24(i+1)/norm(y24(i+1))-y24(i+1);
w24 = w24 + mu24 *e24(i+1)*(signal_ns(:,i+1))';
end;
for k = 0:NumofAntenna-1,
response(k+1,:) = exp(j*k*Kd*sin(theta));
end;
% CALCULATE ARRAY RESPONSE
R = w24*response;
hold on
figure,
plot((theta*180/pi), 20*log10(abs(R)));
title('Beam Pattern');
ylabel('Normalized Beam Pattern(dB) ');
xlabel('Angle(Degrees)');
axis([-90,+90,-100,10]);