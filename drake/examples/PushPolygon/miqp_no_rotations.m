function result = miqp_no_rotations(A_polygon,b_polygon,m,mu,h,r0,p0,rF,N)

checkDependency('gurobi');

A_bar_polygon = [A_polygon(:,2),-A_polygon(:,1)];
num_faces = size(A_polygon,1);

% Decision variables: $$r[n],\dot{r}[n],p[n],\dot{p}[n],\beta_1[n],\beta_2[n],z[n],$$
num_vars_per_timestep = 8 + 3*num_faces;
num_vars = N*num_vars_per_timestep + 6; % r,rdot, and p are also defined for time N+1

% setup helper variable indices 
% ind(n,:) are the indices for the variable at time n
function inds = allTimeSteps(inds_one_timestep,N)
  inds = repmat(inds_one_timestep,N,1) + repmat((0:N-1)'*num_vars_per_timestep,1,size(inds_one_timestep,2));
end
r_inds=allTimeSteps(1:2,N+1); 
rdot_inds=allTimeSteps(3:4,N+1);
p_inds=allTimeSteps(5:6,N+1);
pdot_inds=allTimeSteps(7:8,N);
beta1_inds=allTimeSteps(pdot_inds(1,end)+(1:num_faces),N);
beta2_inds=allTimeSteps(beta1_inds(1,end)+(1:num_faces),N);
z_inds=allTimeSteps(beta1_inds(1,end)+(1:num_faces),N);

model.vtype = [repmat([repmat('C',num_vars_per_timestep-num_faces,1);repmat('B',num_faces,1)],N,1);repmat('C',6,1)];
model.A = sparse(0,num_vars);
model.rhs = [];
model.sense = '';
model.lb = -inf(num_vars,1);
model.ub = inf(num_vars,1);

model.Q = sparse(pdot_inds(:),pdot_inds(:),1+0*pdot_inds(:),num_vars,num_vars);
model.obj = zeros(num_vars,1);

bigM = 1e6;

for n=1:N
  % \sum_i z_i[n] <= 1
  model.A(end+1,z_inds(n,:)) = ones(1,num_faces);
  model.rhs(end+1) = 1;
  model.sense(end+1) = '<';
  
  % m(rdot[n+1]-rdot[n])/h = (A+mu*A_bar)^T\beta1[n] + (A-mu*A_bar)^T beta2[n] 
  cind = size(model.A,1)+(1:2);
  model.A(cind,rdot_inds(n+1,:)) = m/h;
  model.A(cind,rdot_inds(n,:)) = -m/h;
  model.A(cind,beta1_inds(n,:)) = -(A_polygon+mu*A_bar_polygon)';
  model.A(cind,beta2_inds(n,:)) = -(A_polygon-mu*A_bar_polygon)';
  model.rhs(cind) = 0;
  model.sense(cind) = '=';
  
  % r[n+1] = r[n]+h*rdot[n]
  cind = size(model.A,1)+(1:2);
  model.A(cind,r_inds(n+1,:)) = 1;
  model.A(cind,r_inds(n,:)) = -1;
  model.A(cind,rdot_inds(n,:)) = -h;
  model.rhs(cind) = 0;
  model.sense(cind) = '=';
  
  % p[n+1] = p[n]+h*pdot[n-1] 
  cind = size(model.A,1)+(1:2);
  model.A(cind,p_inds(n+1,:)) = 1;
  model.A(cind,p_inds(n,:)) = -1;
  model.A(cind,pdot_inds(n,:)) = -h;
  model.rhs(cind) = 0;
  model.sense(cind) = '=';
  
  % b_polygon-M*(1-z[n]) <= A_polygon p[n]
  cind = size(model.A,1)+(1:num_faces);
  model.A(cind,p_inds(n,:)) = A_polygon;
  model.A(cind,z_inds(n,:)) = -bigM;
  model.rhs(cind) = b_polygon - bigM;
  model.sense(cind) = '>';
  
  % A_polygon p[n] <= b_polygon
  cind = size(model.A,1)+(1:num_faces);
  model.A(cind,p_inds(n,:)) = A_polygon;
  model.rhs(cind) = b_polygon;
  model.sense(cind) = '<';
  
  % 0 \le \beta_1[n],\beta_2[n] \le Mz[n]
  model.lb(beta1_inds(n,:),1) = 0;
  model.lb(beta2_inds(n,:),1) = 0;
  cind = size(model.A,1)+(1:num_faces);
  model.A(cind,beta1_inds(n,:)) = 1;
  model.A(cind,z_inds(n,:)) = -bigM;
  model.rhs(cind) = 0;
  model.sense(cind) = '<';
  cind = size(model.A,1)+(1:num_faces);
  model.A(cind,beta2_inds(n,:)) = 1;
  model.A(cind,z_inds(n,:)) = -bigM;
  model.rhs(cind) = 0;
  model.sense(cind) = '<';
  
end

% r[0] = r0
cind = size(model.A,1)+(1:2);
model.A(cind,r_inds(1,:)) = 1;
model.rhs(cind) = r0;
model.sense(cind) = '=';

% p[0] = p0
cind = size(model.A,1)+(1:2);
model.A(cind,p_inds(1,:)) = 1;
model.rhs(cind) = p0;
model.sense(cind) = '=';

% rdot[0] = 0
cind = size(model.A,1)+(1:2);
model.A(cind,rdot_inds(1,:)) = 1;
model.rhs(cind) = 0;
model.sense(cind) = '=';

% pdot[0] = 0
cind = size(model.A,1)+(1:2);
model.A(cind,pdot_inds(1,:)) = 1;
model.rhs(cind) = 0;
model.sense(cind) = '=';

% r[N+1] = rF
cind = size(model.A,1)+(1:2);
model.A(cind,r_inds(N+1,:)) = 1;
model.rhs(cind) = rF;
model.sense(cind) = '=';

result = gurobi(model);

end
