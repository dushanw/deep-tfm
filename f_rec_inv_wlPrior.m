% 20181107 by Dushan N. Wadduwage

function Xhat = f_rec_inv_wlPrior(pram,Ex,Yhat,gamma,wname)

  % parameters
  if isempty(gamma) 
    gamma = 5e-5;       % the reconstruction is sensitive to the weight of the regularizer (for 64x64 gamma=3e-3 seems to work)
  end
  if isempty(wname) 
    wname       = 'db4';
  end
  
  % make double for sparse;
  Ex    = double(Ex);
  Yhat  = double(Yhat);

  % make A
  i_vec = [1:pram.Ny*pram.Nx*pram.Nt]';
  j_vec = repmat(1:pram.Ny*pram.Nx,[1 pram.Nt])';
  s_vec = Ex(:);
  A = sparse(i_vec,j_vec,s_vec,pram.Ny*pram.Nx*pram.Nt,pram.Ny*pram.Nx);

%  y(find(y(:)<0))=0;  % no negative measurements 
  y = Yhat(:);
  y = y./max(y);

  y = A'*y;           % convert to a (Ny*Nx) by (Ny*Nx) system 
  A = A'*A;

  c = sum(A,1);
  A = A./max(c);
  c = c./max(c);

%    w = diag(sqrt(1./(y+1)));
  w = (sqrt(1./(y+1)));

  L1sum = norm(y(:),1);
  L2sum = sqrt(L1sum); % Target of least square estimation based on Poisson statistics
  eps   = L2sum*1.5;

  %% optimization over x      
%   Psy   = inv(getWaveletmatrices(pram.Ny,pram.Nx,wname));
%   n     = size(A,2);
%   cvx_begin quiet
%       variable x(n)
%       minimize(norm(w.*(A*x-y),2) + gamma*norm(Psy*x,1))        
% 
%       subject to
%           x >= 0;
%           %norm( A * x - y, 2 )<=eps;
%   cvx_end

  %% optimization over alpha  
  Psy   = getWaveletmatrices(pram.Ny,pram.Nx,wname);
  n     = size(Psy,2);
  
  AxPsy = A*Psy;
  cvx_begin 
      variable alph(n)
      % minimize(norm(w.*(A*Psy*alph-y),2) + gamma*norm(alph,1))
      minimize(norm(w.*(AxPsy*alph-y),2) + gamma*norm(alph,1))

      subject to
          Psy*alph >= 0;
          %norm( A * x - y, 2 )<=eps;
  cvx_end
  
  x   = Psy*alph;
  
  %%
%   x(find(isnan(x)))=0;
%   x(find(x(:)<0))=0;

  Xhat = reshape(x,pram.Ny,pram.Nx);
  
  % show reconstruction
  imagesc(imtile([rescale(Xhat) ...
                  rescale(imresize(Yhat(:,:,1,:),[pram.Ny pram.Nx],'nearest'))...
                ]));
  axis image
  
end


function A_waverec2 = getWaveletmatrices(h,w,wname)
  X           = zeros(h,w);
  L           = wmaxlev(size(X),wname);
  [c s]       = wavedec2(X,2,wname);
  A_waverec2  = sparse(h*w,length(c));

  parfor i=1:length(c)
    c1{i}    = zeros(size(c));
    c1{i}(i) = 1;
    A_waverec2(:,i)=sparse(reshape(waverec2(c1{i},s,wname),[h*w,1])); 
  end
  done=1;
end


