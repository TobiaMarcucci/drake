function testManipulatorDynamics

testBrickQuaternion();
testActuatedPendulum();
regressionTestAtlasRPY();
testAtlasQuat();

checkGradients(createFallingBrick('quat'));
checkGradients(createAtlas('rpy'));
checkGradients(createAtlas('quat'));

checkHdotMinus2CoriolisMatrixSkewSymmetricMatrix(createAtlas('rpy'))
end

function robot = createFallingBrick(floating_type)
options.floating = floating_type;
robot = RigidBodyManipulator('FallingBrick.urdf',options);
end

function testBrickQuaternion()
options.floating = 'quat';
m = RigidBodyManipulator('FallingBrick.urdf',options);
nq = m.getNumPositions();
nv = m.getNumVelocities();
q = randn(nq, 1);
v = randn(nv, 1);
[H,C,B] = manipulatorDynamics(m,q,v,false);

valuecheck(H, m.body(2).I, 1e-12);
%TODO: check C
valuecheck([nv 0], size(B));
end

function testActuatedPendulum()
m = RigidBodyManipulator('ActuatedPendulum.urdf');
nq = m.getNumPositions();
nv = m.getNumVelocities();
q = randn(nq, 1);
v = randn(nv, 1);
[H,C,B] = manipulatorDynamics(m,q,v,false);

body = m.body(2);
axis = body.joint_axis;
I = body.I;
H_expected = axis' * I(1:3, 1:3) * axis;
valuecheck(H_expected, H, 1e-12);
% TODO: check C
valuecheck(1, B);
end

function regressionTestAtlasRPY()
replaceMatFile = false;

rng(23415, 'twister');

r = createAtlas('rpy');
nq = r.getNumPositions();
nv = r.getNumVelocities();

nTests = 5;
Hs = cell(nTests, 1);
Cs = cell(nTests, 1);

for i = 1 : nTests
  q = randn(nq, 1);
  v = randn(nv, 1);
  [H, C, B] = manipulatorDynamics(r, q, v, false);
  Hs{i} = H;
  Cs{i} = C;
end

filename = 'regressionTestAtlasManipulatorDynamics.mat';
if replaceMatFile
  save(filename, varname(Hs), varname(Cs), varname(B));
else
  data = load(filename);
  for i = 1 : nTests
    valuecheck(data.Hs{i}, Hs{i}, 1e-10);
    valuecheck(data.Cs{i}, Cs{i}, 1e-10);
  end
  valuecheck(data.B, B);
end

end

function testAtlasQuat()
r = createAtlas('quat');
nv = r.getNumVelocities();

nTests = 5;
for i = 1 : nTests
  q = getRandomConfiguration(r);
  v = randn(nv, 1);
  [H, C, B] = manipulatorDynamics(r, q, v, false);
  kinetic_energy = 1/2 * v' * H * v;
  
  kinsol = r.doKinematics(q, false, false, v);
  kinetic_energy_via_kinsol = computeKineticEnergy(r, kinsol);
  valuecheck(kinetic_energy_via_kinsol, kinetic_energy, 1e-10);
end

end

function out = varname(~)
out = inputname(1);
end

function ret = computeKineticEnergy(manipulator, kinsol)
NB = manipulator.getNumBodies();
ret = 0;

for i = 2 : NB
  twistInBody = relativeTwist(manipulator, kinsol, 1, i, i);
  I = manipulator.body(i).I;
  ret = ret + 1/2 * twistInBody' * I * twistInBody;
end
end

function checkGradients(robot)
q = getRandomConfiguration(robot);
v = randn(robot.getNumVelocities(), 1);
[~,~,~,dH,dC,dB] = manipulatorDynamics(robot, q, v, false);

option.grad_method = 'taylorvar';
[~, ~, ~,dH_geval, dC_geval, dB_geval] = geval(3, @(q, v) manipulatorDynamics(robot, q, v, false), q, v, option);

valuecheck(dH_geval, dH, 1e-10);
valuecheck(dC_geval, dC, 1e-10);
% valuecheck(dB_geval, dB, 1e-10);
end

function checkHdotMinus2CoriolisMatrixSkewSymmetricMatrix(robot)
% Checks that \dot{H} - 2Q(q, v) is skew symmetric, where Q(q, v) is the
% Coriolis matrix. See Lemma 4.2 in Murray, Li, Sastry - A Mathematical
% Introduction to Robotic Manipulation.
% Note that CBar(q, v) = C(q, v) - friction(v) is quadratic in v, i.e. we
% can write the ith entry of Cbar as
% Cbar_i(q, v) = v' * A_i(q) v + N_i(q)
% so dCbar_i(q, v)/dv = 2 * v' * A_i(q)
% We can also write CBar(q, v) as
% CBar_i(q, v) = Q_(i,:)(q, v) * v + N_i(q)
% which shows that Q_(i,:)(q, v) = 1/2 * dCbar_i(q, v)/dv

nq = robot.getNumPositions();
nv = robot.getNumVelocities();
q = getRandomConfiguration(robot);
v = randn(nv, 1);

vToqdot = robot.vToqdot(q);
if nq ~= nv || any(any(vToqdot - eye(nv)))
  error('Hdot - 2Q is only skew symmetric if v = qdot')
end

[~,~,~,dH,dC] = manipulatorDynamics(robot, q, v, false);
qdot = vToqdot * v;
dHdq = dH(:, 1:nq);
Hdot = reshape(dHdq * qdot, [nv, nv]);
[~, dfrictiondv] = robot.computeFrictionForce(v);
dCdv = dC(:, nq + (1:nv)) - dfrictiondv;
coriolis_matrix = 1/2 * dCdv;

Hdot_minus_2_coriolis_matrix = Hdot - 2 * coriolis_matrix;
valuecheck(zeros(nv, nv), Hdot_minus_2_coriolis_matrix + Hdot_minus_2_coriolis_matrix', 1e-10);
end
