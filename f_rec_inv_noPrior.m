% 20200922 by Dushan N. Wadduwage
% Reconstruction by solving the liner system with a psudo inverse
% Yhat = AX => Xhat = A\Yhat

function Xhat = f_rec_inv_noPrior(pram,Ex,Yhat,X0)
  
  % make double for sparse;
  Ex    = double(Ex);
  Yhat  = double(Yhat);

  % make A
  i_vec = [1:pram.Ny*pram.Nx*pram.Nt]';
  j_vec = repmat(1:pram.Ny*pram.Nx,[1 pram.Nt])';
  s_vec = Ex(:);
  A = sparse(i_vec,j_vec,s_vec,pram.Ny*pram.Nx*pram.Nt,pram.Ny*pram.Nx);

  % solve for X using left devide (i.e. psudo inverse)
  
  Xhat = A\Yhat(:);                 % solve original system  
  Xhat = (A'*A)\(A'*Yhat(:));       % solve reduced systems by At
  
  Xhat = reshape(Xhat,pram.Ny,pram.Nx);

  if ~isempty(X0)
    imagesc(imtile([rescale(X0) ...
                    rescale(Xhat) ...
                    rescale(imresize(Yhat(:,:,1,:),[pram.Ny pram.Nx],'nearest'))...
                    ]));
    axis image       
  else
    imagesc(imtile([rescale(Xhat) ...
                    rescale(imresize(Yhat(:,:,1,:),[pram.Ny pram.Nx],'nearest'))...
                    rescale(imresize(mean(Yhat-min(Yhat(:)),3) ,[pram.Ny pram.Nx],'nearest'))...
                  ]));                
    title('Xhat(L),representative Yhat(M),mean Yhat(R)')
    axis image
  end

end