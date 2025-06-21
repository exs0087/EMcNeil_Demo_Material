function dX = odefcn_hw_Q2(t,X,D,S)

    m = 45.93/1000; %kg mass of golf ball
    w = [250 0 0]'; %rad/s.
    g = 9.81; %m/s^2 gravity
    
    dX = zeros(6,1);
  
    dX(1) = X(4); 
    dX(2) = X(5);
    dX(3) = X(6);
    dX(4) = (-D/m)*(X(4)^2) + (S/m)*w(2)*X(6) - (S/m)*w(3)*X(5);
    dX(5) = (-D/m)*(X(5)^2) + (S/m)*w(3)*X(4) - (S/m)*w(1)*X(6);
    dX(6) = (-D/m)*(X(6)^2) + (S/m)*w(1)*X(5) - (S/m)*w(2)*X(4) - g;

end


