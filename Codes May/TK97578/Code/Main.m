%%loading the dataset of weather conditions
% load('DelhiWeatherdata.mat')
% load('DelhiWeatherdatalabels.mat')
load('DelhiWeatherdatan.mat')
load('DelhiWeatherdatatst.mat')
load('classes.mat')
load('classestst.mat')

% %Converting the data from Tabular Form to Matrix Form
% DelhiWeatherdatan = table2array(DelhiWeatherdata);
% DelhiWeatherdatan(1 = DelhiWeatherdatan;

%Labels, 0 - Broken Clouds, 1 - Scatter Clouds, 2 - Light Rain,
%3 - Clear Sky, 4 - Overcast Clouds
%%Separating Features and Labels from the data
x = DelhiWeatherdatan;   %Features
y = classes;            %Labels

%%Classification using Multi-Class SVM
%Training
Mdl = fitcecoc(x,y);

%Testing
inp = input('Enter any number between 1 and 9: ');
Xnew = DelhiWeatherdatatst(inp,:);
szdwdt = size(DelhiWeatherdatatst);
ypred = predict(Mdl,Xnew);  %prediction

%Predicted output
if ypred{1} == '0'
    out = 'Broken Clouds'
elseif ypred{1} == '1'
    out = 'Scattered Clouds'
elseif ypred{1} == '2'
    out = 'Light Rain'
elseif ypred{1} == '3'
    out = 'Clear Sky'
elseif ypred{1} == '4'
    out = 'Overcast Clouds'
end

%accuracy
cnt = 0;
for in = 1:szdwdt(1)
    ypredn = predict(Mdl,DelhiWeatherdatatst(in,:));
    if ypredn{1} == classestst{in}
        cnt = cnt+1;
    end
end

accuracy = (cnt/szdwdt(1))*100;
