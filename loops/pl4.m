%Switch case in the matlab
a=input('enter the day in the week:','s')
switch a
    case'monday'
        disp('Day 1')
    case'tuesday'
        disp('Day 2')
    case'wednesday'
        disp('Day 3')
    case'thursday'
        disp('Day 4')
    case'firday'
        disp('Day 5')
    case'saturday'
        disp('weekend')
    case'sunday'
        disp('holiday')
    otherwise
        disp('enter the valid day')
end
