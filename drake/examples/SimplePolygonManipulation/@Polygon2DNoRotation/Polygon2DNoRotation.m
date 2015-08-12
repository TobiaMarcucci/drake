classdef Polygon2DNoRotation < DrakeSystem
  
  properties
    A,b;      % convex polygon is defined as Ax \le b
    m = 1;
    mu = .1;  % coefficient of friction
    h = .1;   % timestep
  end

  methods
    function obj = Polygon2DNoRotation(A,b)
      obj = obj@DrakeSystem(0,6,2,4,false,true);
      obj.A = A;
      obj.b = b;
    end
    
    function xn = update(obj,t,x,u)
      error('todo: implement this');
    end
    
    function planDemo(obj)
      r0=zeros(2,1); 
      p0=zeros(2,1);
      rF=randn(2,1);
      N=10;
      v = Polygon2DNoRotationVisualizer(obj);
      v.drawWrapper(0,[r0;p0]);
      
      ytraj=miqpPlanner(obj,r0,p0,rF,N);
      v.playback(ytraj);
    end
  end
  
  methods(Static)
    
    function obj = box()
      A = [1,0;0,1;-1,0;0,-1]; b=[1;1;1;1];
      obj = Polygon2DNoRotation(A,b);
    end
    
    function obj = random(num_vertices)
      vertices = randn(num_vertices,2);
      vertices = vertices - repmat(mean(vertices,1),size(vertices,1),1);
      [A,b] = vert2lcon(vertices);
      obj = Polygon2DNoRotation(A,b);
    end
  end
end