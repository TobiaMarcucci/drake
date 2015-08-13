classdef Polygon2DNoRotationVisualizer < Visualizer

  properties
    vertices
  end

  methods
    function obj = Polygon2DNoRotationVisualizer(push_poly)
      typecheck(push_poly,'Polygon2DNoRotation');
      obj = obj@Visualizer(push_poly.getOutputFrame());
      obj.vertices = lcon2vert(push_poly.A,push_poly.b)';
      k = convhull(obj.vertices(1,:),obj.vertices(2,:));
      obj.vertices = obj.vertices(:,k);
    end
    
    function draw(obj,~,y)
      r = y(1:2);
      p = y(3:4);
      active_face = y(5);
      pts = obj.vertices + repmat(r,1,size(obj.vertices,2));
      patch(pts(1,:),pts(2,:),'c');
      if (active_face>.01)
        plot(r(1)+p(1),r(2)+p(2),'g.','MarkerSize',12);
      else
        plot(r(1)+p(1),r(2)+p(2),'r.','MarkerSize',12);
      end
      axis equal;
      if isempty(obj.axis)
        axis([-4 4 -4 4]);
      else
        axis(obj.axis);
      end
    end
  end
end