classdef Polygon2D < DrakeSystem
  
  properties
    A,b;      % convex polygon is defined as Ax \le b
    m = 1;
    mu = .1;  % coefficient of friction
    h = .1;   % timestep
    c = .2;   % damping coefficient
  end

  methods
    function obj = Polygon2D(A,b)
      obj = obj@DrakeSystem(0,8,2,6,false,true);
      obj = obj.setOutputFrame(Polygon2D.singletonOutputFrame());
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
      v.drawWrapper(0,[r0;p0;0]);
      v.playback_speed = 0.25;
      
      ytraj=miqpPlanner(obj,r0,p0,rF,N);
      v.playback(ytraj);
    end
  end
  
  methods(Static)
    
    function obj = box()
      A = [1,0;0,1;-1,0;0,-1]; b=[1;1;1;1];
      obj = Polygon2D(A,b);
    end
    
    function obj = random(num_vertices)
      vertices = randn(num_vertices,2);
      vertices = vertices - repmat(mean(vertices,1),size(vertices,1),1);
      [A,b] = vert2lcon(vertices);
      obj = Polygon2D(A,b);
    end
    
    function fr = singletonOutputFrame()
      fr = SingletonCoordinateFrame('PolygonPosition',6,'x',{'rx','ry','px','py','theta','active_face'});
    end
  end
end