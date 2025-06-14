function [y1,y2,y3]=sat(s)
%Erin McNeil

eps1 = 0;
if s(1) > eps1
    y1=1;
elseif norm(s(1))<= eps1
    y1=s(1)/eps1;
elseif s(1) < -eps1
    y1=-1;
else
    errormsg("something's amiss")
end

eps2 = 0;
if s(2) > eps2
    y2=1;
elseif norm(s(2))<= eps2
    y2=s(2)/eps2;
elseif s(2) < -eps2
    y2=-1;
else
    errormsg("something's amiss")
end

eps3 = 0;
if s(3) > eps3
    y3=1;
elseif norm(s(3))<= eps3
    y3=s(3)/eps3;
elseif s(3) < -eps3
    y3=-1;
else
    errormsg("something's amiss")
end