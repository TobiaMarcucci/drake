function ytraj = imiqpPlanner(obj,r0,p0,theta0,rF,thetaF,N)

checkDependency('gurobi');
%checkDependency('spotless');  % only used below for debugging constraints

A_bar = [obj.A(:,2),-obj.A(:,1)];
num_faces = size(obj.A,1);

%% Decision variables: $$r[n],\dot{r}[n],p[n],\dot{p}[n],\beta_1[n],\beta_2[n],z[n],$$
num_vars_per_timestep = 10 + 3*num_faces;
num_vars = N*num_vars_per_timestep + 8; % r,rdot, and p,theta,thetadot are also defined for time N+1

%% setup helper variable indices 
% ind(n,:) are the indices for the variable at time n
%x = msspoly('x',num_vars);
function inds = allTimeSteps(inds_one_timestep,N,name)
  inds = repmat(inds_one_timestep,N,1) + repmat((0:N-1)'*num_vars_per_timestep,1,size(inds_one_timestep,2));
%  if nargin>2
%    for i=1:N
%      x(inds(i,:)) = msspoly([name,char('a'+(i-1))],length(inds_one_timestep));
%    end
%  end
end
r_inds=allTimeSteps(1:2,N+1,'r'); 
p_inds=allTimeSteps(3:4,N+1,'p');
theta_inds=allTimeSteps(5,N+1,'th');
rdot_inds=allTimeSteps(6:7,N+1,'dr');
thetadot_inds=allTimeSteps(8,N+1,'dth');
pdot_inds=allTimeSteps(9:10,N,'dp');
beta1_inds=allTimeSteps(10+(1:num_faces),N,'ba');
beta2_inds=allTimeSteps(beta1_inds(1,end)+(1:num_faces),N,'bb');
z_inds=allTimeSteps(beta2_inds(1,end)+(1:num_faces),N,'z');

%% setup model
model.vtype = [repmat([repmat('C',num_vars_per_timestep-num_faces,1);repmat('B',num_faces,1)],N,1);repmat('C',8,1)];
model.A = sparse(0,num_vars);
model.rhs = [];
model.sense = '';
model.lb = -inf(num_vars,1);
model.ub = inf(num_vars,1);

model.obj = zeros(num_vars,1);

bigM = 1e2;

model.Q = sparse(num_vars,num_vars);
result.x = zeros(num_vars,1);
x0 = ones(num_vars,1);  % initialize iterations

fn_beta1 = -(obj.A+obj.mu*A_bar)';
fn_beta2 = -(obj.A-obj.mu*A_bar)';

iter=0;
while norm(x0 - result.x,2) > 1e-3
  iter=iter+1;
  disp(['iter ',num2str(iter)]);drawnow;
  x0 = result.x;

  for n=1:N
    R0 = rotmat(x0(theta_inds(n,:)));
    dR0 = drotmat(x0(theta_inds(n,:)));
    
    model.Q(pdot_inds(n,:),pdot_inds(n,:)) = R0'*R0;
    model.Q(rdot_inds(n,:),rdot_inds(n,:)) = eye(2);
    model.Q(thetadot_inds(n,:),thetadot_inds(n,:)) = norm(dR0*x0(pdot_inds(n,:)));
    
    %% sum_i z_i[n] <= 1
    model.A(end+1,z_inds(n,:)) = ones(1,num_faces);
    model.rhs(end+1) = 1;
    model.sense(end+1) = '<';
    
    f0n = fn_beta1*x0(beta1_inds(n,:)) + fn_beta2*x0(beta2_inds(n,:));
    
    %% m\frac{\dot{r}[n+1]-\dot{r}[n]}{h} =& R(\theta_0[n])f[n] -\frac{\partial R(\theta_0[n])}
    cind = size(model.A,1)+(1:2);
    model.A(cind,rdot_inds(n+1,:)) = obj.m/obj.h*eye(2);
    model.A(cind,rdot_inds(n,:)) = (-obj.m/obj.h + obj.c)*eye(2);
    model.A(cind,beta1_inds(n,:)) = -R0*fn_beta1;
    model.A(cind,beta2_inds(n,:)) = -R0*fn_beta2;
    model.A(cind,theta_inds(n,:)) = -dR0*f0n;
    model.rhs(cind) = 0;
    model.sense(cind) = '=';
      
    %% I\frac{\dot\theta[n+1]-\dot\theta[n]}{h} =& p_0[n] \times f[n] + p[n] \times f_0[n] - c \dot\theta[n]
    cind = size(model.A,1)+1;
    model.A(cind,thetadot_inds(n+1,:)) = obj.I/obj.h;
    model.A(cind,thetadot_inds(n,:)) = -obj.I/obj.h + obj.c;
    model.A(cind,beta1_inds(n,:)) = -x0(p_inds(n,1))*fn_beta1(2,:) + x0(p_inds(n,2))*fn_beta1(1,:);
    model.A(cind,beta2_inds(n,:)) = -x0(p_inds(n,1))*fn_beta2(2,:) + x0(p_inds(n,2))*fn_beta2(1,:);
    model.A(cind,p_inds(n,1)) = f0n(2);
    model.A(cind,p_inds(n,2)) = -f0n(1);
    model.rhs(cind) = 0;
    model.sense(cind) = '=';

    for i=1:2
      %% r[n+1] = r[n]+h*rdot[n]
      cind = size(model.A,1)+1;
      model.A(cind,r_inds(n+1,i)) = 1;
      model.A(cind,r_inds(n,i)) = -1;
      model.A(cind,rdot_inds(n,i)) = -obj.h;
      model.rhs(cind) = 0;
      model.sense(cind) = '=';
      
      %% p[n+1] = p[n]+h*pdot[n]
      cind = size(model.A,1)+1;
      model.A(cind,p_inds(n+1,i)) = 1;
      model.A(cind,p_inds(n,i)) = -1;
      model.A(cind,pdot_inds(n,i)) = -obj.h;
      model.rhs(cind) = 0;
      model.sense(cind) = '=';
    end
    
    %% theta[n+1] = theta[n]+h*thetadot[n]
    cind = size(model.A,1)+1;
    model.A(cind,theta_inds(n+1,:)) = 1;
    model.A(cind,theta_inds(n,:)) = -1;
    model.A(cind,thetadot_inds(n,:)) = -obj.h;
    model.rhs(cind) = 0;
    model.sense(cind) = '=';
    
    %% b_polygon-M*(1-z[n]) <= obj.A p[n]
    cind = size(model.A,1)+(1:num_faces);
    model.A(cind,p_inds(n,:)) = -obj.A;
    model.A(cind,z_inds(n,:)) = bigM*eye(num_faces);
    model.rhs(cind) = -obj.b + bigM;
    model.sense(cind) = '<';
    
    %% obj.A p[n] <= b_polygon
    cind = size(model.A,1)+(1:num_faces);
    model.A(cind,p_inds(n,:)) = obj.A;
    model.rhs(cind) = obj.b;
    model.sense(cind) = '<';
    
    %% 0 \le \beta_1[n],\beta_2[n] \le Mz[n]
    model.lb(beta1_inds(n,:),1) = 0;
    model.lb(beta2_inds(n,:),1) = 0;
    cind = size(model.A,1)+(1:num_faces);
    model.A(cind,beta1_inds(n,:)) = eye(num_faces);
    model.A(cind,z_inds(n,:)) = -bigM*eye(num_faces);
    model.rhs(cind) = 0;
    model.sense(cind) = '<';
    cind = size(model.A,1)+(1:num_faces);
    model.A(cind,beta2_inds(n,:)) = eye(num_faces);
    model.A(cind,z_inds(n,:)) = -bigM*eye(num_faces);
    model.rhs(cind) = 0;
    model.sense(cind) = '<';
  end
  
  %% r[0] = r0
  cind = size(model.A,1)+(1:2);
  model.A(cind,r_inds(1,:)) = eye(2);
  model.rhs(cind) = r0;
  model.sense(cind) = '=';
  
  %% p[0] = p0
  cind = size(model.A,1)+(1:2);
  model.A(cind,p_inds(1,:)) = eye(2);
  model.rhs(cind) = p0;
  model.sense(cind) = '=';
  
  %% theta[0] = theta0
  cind = size(model.A,1)+1;
  model.A(cind,theta_inds(1,:)) = 1;
  model.rhs(cind) = theta0;
  model.sense(cind) = '=';
  
  %% rdot[0] = 0
  cind = size(model.A,1)+(1:2);
  model.A(cind,rdot_inds(1,:)) = eye(2);
  model.rhs(cind) = 0;
  model.sense(cind) = '=';
  
  %% pdot[0] = 0
  cind = size(model.A,1)+(1:2);
  model.A(cind,pdot_inds(1,:)) = eye(2);
  model.rhs(cind) = 0;
  model.sense(cind) = '=';
  
  %% thetadot[0] = 0
  cind = size(model.A,1)+1;
  model.A(cind,thetadot_inds(1,:)) = 1;
  model.rhs(cind) = 0;
  model.sense(cind) = '=';
  
  
  %% r[N+1] = rF
  cind = size(model.A,1)+(1:2);
  model.A(cind,r_inds(N+1,:)) = eye(2);
  model.rhs(cind) = rF;
  model.sense(cind) = '=';
  
  %% theta[N+1] = thetaF
  cind = size(model.A,1)+1;
  model.A(cind,theta_inds(N+1,:)) = 1;
  model.rhs(cind) = thetaF;
  model.sense(cind) = '=';

  %[model.A * x - model.rhs']
  
  %% solve
  result = gurobi(model);
end

r = result.x(r_inds)';
p = result.x(p_inds)';
theta = result.x(theta_inds)';
rdot = result.x(rdot_inds)';
pdot = result.x(pdot_inds)';
thetadot = result.x(thetadot_inds)';
beta1 = result.x(beta1_inds)';
beta2 = result.x(beta2_inds)';
z = result.x(z_inds)';

active_face = (1:num_faces)*z;

ytraj = DTTrajectory((0:N)*obj.h,[r;p;0*p(1,:);[active_face,0]]);
ytraj = setOutputFrame(ytraj,getOutputFrame(obj));

end
