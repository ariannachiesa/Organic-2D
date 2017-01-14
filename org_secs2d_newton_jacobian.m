## Assemble jacobian matrix for Newton's method.
function jacobian = org_secs2d_newton_jacobian ...
                    (device, material, constants,
                     quadrature, charge_n,
                     algorithm, V, n, F,
                     deltat, A, B, C, r, _indexing)
  
  nnodes = columns(device.msh.p);
  assert (numel (unique (device.pins)) == 2);
  nelements = columns(device.msh.t);
  
  numcontacts = numel(device.pins);
  ndofs = 2 * nnodes + numel (F) + numcontacts;
  
  if (! isstruct (_indexing))
    indexing.V = _indexing (1:nnodes);
    indexing.n = _indexing ((1:nnodes) + nnodes);
    indexing.F = _indexing (2*nnodes + (1:numel (F)));
    indexing.I = _indexing ((2*nnodes + numel (F) + 1) : end);
  else
    indexing = _indexing;
  endif
  
  persistent dnodes = device.dnodes;
  persistent alldnodes = device.alldnodes;
  
  intnodes = setdiff (1:nnodes, alldnodes);
  intdofsV = indexing.V (intnodes);
  intdofsn = indexing.n (intnodes);
  
  ## COMPUTING COEFFICIENTS
  epsilon = material.eps_semic * ones (nelements, 1);
  if isfield(device, 'insulator')
    epsilon(device.insulator) = material.eps_ins;
  endif
  
  ## COMPUTING FIRST ROW
  MM = bim2a_reaction (device.msh, ! device.insulator, 1);
  M = bim2a_reaction (device.msh, 1, 1);
  
  A11 = bim2a_laplacian (device.msh, epsilon, 1);
  A12 = +constants.q * MM;
  
  ## Physical models.
  [mobility, alpha, der] = org_physical_models2d ...
                                (n, material, device, constants,
                                 quadrature, charge_n);
  ##
  jacobian = sparse (ndofs, ndofs);
  
  ## ASSEMBLING FIRST ROW
  jacobian(intdofsV, indexing.V) += A11(intnodes, :);
  jacobian(intdofsV, indexing.n) += A12(intnodes, :);
  
  ## ENFORCING BCs ON FIRST ROW
  diagM = diag(M);
  for ii = 1:numcontacts
    ddofsV = indexing.V(dnodes{ii});
    pindof = indexing.F(device.pins(ii));
    
    vec = diagM(dnodes{ii});
    
    jacobian(ddofsV, ddofsV) += diag(vec);
    jacobian(ddofsV, pindof) -= vec;
  endfor

  ## COMPUTING SECOND ROW
  A21 = -bim2a_laplacian (device.msh, mobility.n, n);
  A22 = sparse (nnodes, nnodes);
  
  [nx, ny] = bim2c_pde_gradient(device.msh,n);  
  gradn = [nx; ny];
  
  [Vx, Vy] = bim2c_pde_gradient(device.msh,V);  
  gradV = [Vx; Vy];
  
  fluxn = - (der.dalpha.n' .* gradn - gradV / constants.Vth) ;
  
  A22(device.scnodes, device.scnodes) = ...
  bim2a_advection_diffusion ...
    (device.redmsh, mobility.n(! device.insulator) * constants.Vth,
     1, alpha.n(device.scnodes), fluxn(:, ! device.insulator));
  
  A22 += bim2a_reaction (device.msh, ! device.insulator, 1 / deltat);
  
  ## ASSEMBLING SECOND ROW
  jacobian(intdofsn, indexing.V) += A21(intnodes, :);
  jacobian(intdofsn, indexing.n) += A22(intnodes, :);
  
  ## ENFORCING BCs ON SECOND ROW
  for ii = 1:numcontacts
    ddofsn = indexing.n(dnodes{ii});
        
    vec = diagM(dnodes{ii});
    
    if (ii == 1) # Metal/semic. interface.
        jacobian(ddofsn, ddofsn) += diag (vec);
    endif
  endfor
  
  ## ADJUST FOR ZERO INSULATOR CHARGE
  if isfield (device, 'scnodes')
    insn = !device.scnodes;
    jacobian(indexing.n(insn), :) = 0;
    jacobian(indexing.n(insn), indexing.n(insn)) += diag (diagM(insn));
  endif
  
  ## ASSEMBLING THIRD ROW
  jacobian(indexing.F, indexing.F) = (A / deltat) + B;
  jacobian(indexing.F, indexing.I) = r;
  
  ## COMPUTING FOURTH ROW
  jacobian(indexing.I, indexing.I) = eye (numel (indexing.I));
  
  for ii = 1:numcontacts
    rr = dnodes{device.pins(ii)};
    
    # Displacement current.
    jacobian(indexing.I(ii), indexing.V) -= ...
      device.section * sum (A11(rr,:), 1) / deltat;
    jacobian(indexing.I(ii), indexing.n) -= ...
      device.section * sum (A12(rr,:), 1) / deltat;
    
    # Electron current.
    jacobian(indexing.I(ii), indexing.V) -= ...
      - device.section * constants.q * sum (A21(rr,:), 1);
    jacobian(indexing.I(ii), indexing.n) -= ...
      - device.section * constants.q * sum (A22(rr,:), 1);
  endfor
  
  jacobian(indexing.V, :) /= algorithm.rowscaling(1);
  jacobian(indexing.n, :) /= algorithm.rowscaling(2);
  jacobian(indexing.F, :) /= algorithm.rowscaling(3);
  jacobian(indexing.I, :) /= algorithm.rowscaling(4);
  
  jacobian(:, indexing.V) *= algorithm.colscaling(1);
  jacobian(:, indexing.n) *= algorithm.colscaling(2);
  jacobian(:, indexing.F) *= algorithm.colscaling(3);
  jacobian(:, indexing.I) *= algorithm.colscaling(4);
  
endfunction
