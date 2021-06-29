function dydt = Kuramoto_model(t,y,Nosc,myomega,coupling)
% Ref: http://classes.soe.ucsc.edu/ams214/Winter09/lecturenotes/Kuramoto.m
% This function simulates the Kuramoto model.
dydt = zeros(Nosc,1);

% Instrinsic part
for i = 1 : Nosc
    dydt(i) = myomega(i);
    for j = 1 : Nosc
        dydt(i) = dydt(i) + (coupling/Nosc)*sin(2*pi*(y(j)-y(i)));
    end
    dydt(i) = dydt(i) + randn;
end




