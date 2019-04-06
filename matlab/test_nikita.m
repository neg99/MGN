clear all; close all;
addpath ~/repos/slra;

%Ns = [20];

Ns = round(exp(linspace(log(20),log(50000),20))');

settings = {0, 'p'; ...
            1, 'l';
            1, 'p';
            };
    



dist = zeros(length(Ns), size(settings,1),2);

methodnames = cell(size(settings,1),1);

for j=1:size(settings,1);
  methodnames{j} = sprintf('corr = %d, meth = %s', ...
                   settings{j,1}, settings{j,2});
                 
  disp(methodnames{j});           
    
  % Set common settings          
  opt = struct('ls_correction', settings(j,1), 'method', settings(j,2));
  opt.disp = 'notify';
  opt.epsabs = 0;
  opt.epsrel = 0;
  opt.tol = 0;
  opt.epsgrad = 1e-14;
  opt.epsx = 0;
  opt.maxiter = 300;
  opt.reggamma = 1e-12;
  s.m = 4;
  a0 = [1, -3, 3, -1];
  a = a0 + 1e-6 * ones(1,4);
  opt.Rini = a;

      
  for i=1:length(Ns);
    N = Ns(i)  
    x = linspace(-1,1,N)';

    % Generate the series
    YNhat = x.^2; YNhat = YNhat / norm(YNhat);
    NNhat = abs(x); NNhat = NNhat / norm(NNhat);
    V = zeros(N,6);
    for n=1:size(V,2)
      Z = legendre(n-1,x);
      V(:,n) = Z(1,:)';
    end  
    
    %V = x.^(0:5);
    [Q,R] =qr(V, 0);
    NN = NNhat - Q * (Q' * NNhat);
    %[Q1,R] = qr([V,NNhat], 0);
    %NN = Q1(:,7) * R(7,7);
    XN = YNhat + NN;
  
    [y, info] = slra(XN, struct('m',4), 3, opt);
    dist(i,j,1) = norm(YNhat - y);
    dist(i,j,2) = norm(conv(y, info.Rh, 'valid')) / norm(info.Rh);
  end  
end
h = figure;

subplot(121);
loglog(Ns, dist(:,:,1), 'LineWidth', 2);
xlabel('N');
ylabel('Distance to solution');
subplot(122);
loglog(Ns, dist(:,:,2), 'LineWidth', 2);
ylabel('Discrepancy');
xlabel('N');
legend(methodnames{:}, 'location', 'SouthEast');
export_fig('slra_matlab_results.pdf', '-pdf');
