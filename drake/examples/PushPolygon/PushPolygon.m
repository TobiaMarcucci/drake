classdef PushPolygon < DrakeSystem
  
  properties
    A,b;      % polygon is defined as Ax \le b
    mu = .1;  % coefficient of friction
    h = .1;   % timestep
  end

  methods
    
  end
end