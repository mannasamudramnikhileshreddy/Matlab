answer = inputdlg('Enter a number:',...
             'Sample', [1 50])
a = str2num(answer{1})
if a<0
   f = msgbox(["Number is less than zero"]);
elseif a>0
    f = msgbox(["Number is greater than zero"]);
else
   f = msgbox(["Number is greater than zero"]);
end

