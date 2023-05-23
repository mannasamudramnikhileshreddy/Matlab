function Energy = Channel_allocate_FIL(H,D)
%Number of users
Kr = size(H,1);

%Total number of antennas
N = size(H,2);

%If D matrix is not provided, all antennas can transmit to everyone
if nargin<2
    D = repmat( eye(N), [1 1 Kr]);
end

%Energy Harvest
Energy = zeros(size(H'));

for k = 1:Kr
    effectivechannel = (H*D(:,:,k))'; %Effective channels
    channelinversion = effectivechannel/(effectivechannel'*effectivechannel);
    Energy(:,k) = channelinversion(:,k)/norm(channelinversion(:,k));  %Normalization 
end
