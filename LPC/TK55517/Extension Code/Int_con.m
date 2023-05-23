function Interference = CA_fAll(H,D)

%Number of users
Kr = size(H,1);

%Total number of antennas
N = size(H,2);

%If D matrix is not provided, all antennas can transmit to everyone
if nargin<2
    D = repmat( eye(N), [1 1 Kr]);
end

% location of Device
Interference = zeros(size(H'));

%Interference power constraint
for k = 1:Kr
    channelvector = (H(k,:)*D(:,:,k))'; %Useful channel
    Interference(:,k) = channelvector/norm(channelvector); %Normalization of useful channel devices
end
