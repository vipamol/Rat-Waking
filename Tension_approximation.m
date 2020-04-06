function [Ten] = Tension_approximation(x,xdot,dt,Kpe,Kse,C,T)
if isnan(T(1))
    T(1) = 0;
end

   for i = 1:length(x)
       dT = Kse*Kpe/C*x(i) + Kse*xdot(i) - T(i)*(Kse + Kpe)/C;
       T(i+1) = T(i)+ dT*dt;
       if T(i+1) < 0
           T(i+1) = 0;
       end
   end
      Ten = T;
end
