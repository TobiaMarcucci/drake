function [twist, dtwist] = relativeTwist(obj, kinsol, base, end_effector, expressed_in)
% RELATIVETWIST Computes the relative twist between base and end effector
% @param transforms homogeneous transforms from link to world (usually
% obtained from doKinematics as kinsol.T)
% @param kinsol solution structure returned by doKinematics
% @param twists twists of links with respect to world, expressed in world
% @param base index of rigid body that will be considered the base
% @param end_effector index of rigid body that will be considered the end
% effector
% @param expressed_in index of rigid body in whose frame the end result
% will be expressed
% @retval relative twist of end_effector with respect to base, expressed in
% expressed_in

compute_gradient = nargout > 1;

twist = kinsol.twists{end_effector} - kinsol.twists{base};
if compute_gradient
  dtwistdq = kinsol.dtwistsdq{end_effector} - kinsol.dtwistsdq{base};
  
  dtwistdv = zeros(numel(twist), length(kinsol.v)) * kinsol.v(1);
  [J, v_indices] = geometricJacobian(obj, kinsol, base, end_effector, expressed_in);
  dtwistdv(:, v_indices) = J;
end

if expressed_in ~= 1
  T = kinsol.T{expressed_in};
  Tinv = homogTransInv(T);
  if compute_gradient
    dT = kinsol.dTdq{expressed_in};
    dTinv = dinvT(T, dT);
    dtwistdq = dAdHTimesX(Tinv, twist, dTinv, dtwistdq);
  end
  twist = transformTwists(Tinv, twist);
end

if compute_gradient
  dtwist = [dtwistdq dtwistdv];
end

end