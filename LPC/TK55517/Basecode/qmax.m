function Qmax = Channel_allocate_PIL(H,eta,D)

%Number of users
Kr = size(H,1);

%Total number of antennas
N = size(H,2);

%If eta vector is not provided, all values are set to unity
if nargin<2
    eta = ones(Kr,1);
end

%If D matrix is not provided, all antennas can transmit to everyone
if nargin<3
    D = repmat( eye(N), [1 1 Kr]);
end

%location of devices

Qmax = zeros(size(H'));


for k = 1:Kr
    effectivechannel = (H*D(:,:,k))'; %Effective channels
    projectedchannel = (eye(N)/eta(k)+effectivechannel*effectivechannel')\effectivechannel(:,k); %Compute allocation of channel inversion
    Qmax(:,k) = projectedchannel/norm(projectedchannel);  %Normalization 
end
