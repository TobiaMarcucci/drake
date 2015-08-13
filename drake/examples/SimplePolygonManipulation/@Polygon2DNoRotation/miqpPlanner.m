function ytraj = miqpNoRotations(obj,r0,p0,rF,N)

checkDependency('gurobi');

A_bar = [obj.A(:,2),-obj.A(:,1)];
num_faces = size(obj.A,1);

%% Decision variables: $$r[n],\dot{r}[n],p[n],\dot{p}[n],\beta_1[n],\beta_2[n],z[n],$$
num_vars_per_timestep = 8 + 3*num_faces;
num_vars = N*num_vars_per_timestep + 6; % r,rdot, and p are also defined for time N+1

%% setup helper variable indices 
% ind(n,:) are the indices for the variable at time n
x = msspoly('x',num_vars);
function inds = allTimeSteps(inds_one_timestep,N,name)
  inds = repmat(inds_one_timestep,N,1) + repmat((0:N-1)'*num_vars_per_timestep,1,size(inds_one_timestep,2));
  if nargin>2
    for i=1:N
      x(inds(i,:)) = msspoly([name,char('a'+(i-1))],length(inds_one_timestep));
    end
  end
end
r_inds=allTimeSteps(1:2,N+1,'r'); 
rdot_inds=allTimeSteps(3:4,N+1,'rd');
p_inds=allTimeSteps(5:6,N+1,'p');
pdot_inds=allTimeSteps(7:8,N,'pd');
beta1_inds=allTimeSteps(pdot_inds(1,end)+(1:num_faces),N,'ba');
beta2_inds=allTimeSteps(beta1_inds(1,end)+(1:num_faces),N,'bb');
z_inds=allTimeSteps(beta2_inds(1,end)+(1:num_faces),N,'z');

%% setup model
model.vtype = [repmat([repmat('C',num_vars_per_timestep-num_faces,1);repmat('B',num_faces,1)],N,1);repmat('C',6,1)];
model.A = sparse(0,num_vars);
model.rhs = [];
model.sense = '';
model.lb = -inf(num_vars,1);
model.ub = inf(num_vars,1);

cost_inds = [pdot_inds(:);rdot_inds(:)];
model.Q = sparse(cost_inds,cost_inds,1+0*cost_inds,num_vars,num_vars);
model.obj = zeros(num_vars,1);

bigM = 1e2;

for n=1:N
  %% sum_i z_i[n] <= 1
  model.A(end+1,z_inds(n,:)) = ones(1,num_faces);
  model.rhs(end+1) = 1;
  model.sense(end+1) = '<';
  
  for i=1:2
    %% m(rdot[n+1]-rdot[n])/h = -(A+mu*A_bar)^T\beta1[n] + -(A-mu*A_bar)^T beta2[n] - c rdot[n]
    cind = size(model.A,1)+1;
    model.A(cind,rdot_inds(n+1,i)) = obj.m/obj.h;
    model.A(cind,rdot_inds(n,i)) = -obj.m/obj.h + obj.c;
    model.A(cind,beta1_inds(n,:)) = (obj.A(:,i)+obj.mu*A_bar(:,i))';
    model.A(cind,beta2_inds(n,:)) = (obj.A(:,i)-obj.mu*A_bar(:,i))';
    model.rhs(cind) = 0;
    model.sense(cind) = '=';
  
    %% r[n+1] = r[n]+h*rdot[n]
    cind = size(model.A,1)+1;
    model.A(cind,r_inds(n+1,i)) = 1;
    model.A(cind,r_inds(n,i)) = -1;
    model.A(cind,rdot_inds(n,i)) = -obj.h;
    model.rhs(cind) = 0;
    model.sense(cind) = '=';
  
    %% p[n+1] = p[n]+h*pdot[n-1]
    cind = size(model.A,1)+1;
    model.A(cind,p_inds(n+1,i)) = 1;
    model.A(cind,p_inds(n,i)) = -1;
    model.A(cind,pdot_inds(n,i)) = -obj.h;
    model.rhs(cind) = 0;
    model.sense(cind) = '=';
  end
    
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

%% r[N+1] = rF
cind = size(model.A,1)+(1:2);
model.A(cind,r_inds(N+1,:)) = eye(2);
model.rhs(cind) = rF;
model.sense(cind) = '=';

%[model.A * x - model.rhs']

%% solve
result = gurobi(model);

r = result.x(r_inds)';
rdot = result.x(rdot_inds)';
p = result.x(p_inds)';
pdot = result.x(pdot_inds)';
beta1 = result.x(beta1_inds)';
beta2 = result.x(beta2_inds)';
z = result.x(z_inds)';

active_face = (1:num_faces)*z;

ytraj = DTTrajectory((0:N)*obj.h,[r;p;[active_face,0]]);
ytraj = setOutputFrame(ytraj,getOutputFrame(obj));

end
