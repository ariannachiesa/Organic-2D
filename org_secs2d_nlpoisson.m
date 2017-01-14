## Solve the non-linear Poisson problem using Newton's algorithm.
function [phiout, nout, resnrm, iter] = ...
         org_secs2d_nlpoisson (device, material,
                               constants, algorithm,
                               phi0, charge_n)
  
  nnodes = columns(device.msh.p);
  nelements  = columns(device.msh.t);

  intnodes = setdiff (1:nnodes, device.alldnodes);
  
  # Assemble system matrices.
  epsilon = material.eps_semic * ones (nelements, 1);
  if isfield(device, 'insulator')
    epsilon(device.insulator) = material.eps_ins;
  endif

  A = bim2a_laplacian (device.msh, epsilon, 1);
  M = bim2a_reaction (device.msh, ! device.insulator, 1);
  
  # Newton's algorithm.
  phi    = phi0;
  phiout = phi0;

  for iter = 1:algorithm.pmaxit
      phiout = phi;
            
      [rho, drho] = charge_n (phiout);
      
      res = (A * phiout - M * rho);
      
      jac = A - M .* diag (drho);
      
      dphi = -jac(intnodes, intnodes) \ res(intnodes);
      
      resnrm(iter) = norm (dphi, inf);
      
      #dphi (dphi >   algorithm.ptoll) =   algorithm.ptoll;
      #dphi (dphi < - algorithm.ptoll) = - algorithm.ptoll;
      phi(intnodes) += dphi;
      
      if resnrm(iter) < algorithm.ptoll
        break
      endif
  endfor

  phiout = phi;
  
  # Post-processing.
  nout = zeros(nnodes,1);
  nout(device.scnodes) = -charge_n (phiout(device.scnodes)) / constants.q;
endfunction
