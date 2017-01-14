## Assemble residual vector for Newton's method.
function res = org_secs2d_newton_residual ...
               (device, material, constants,
                quadrature, charge_n, algorithm,
                V, n, F, I, V0, n0, F0, I0,
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
  
  intnodes = setdiff (1:nnodes, alldnodes); # internal nodes
  indexV = 1:2:2*nnodes;
  indexn = 2:2:2*nnodes;
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
  [mobility, alpha] = org_physical_models2d ...
                        (n, material, device, constants,
                         quadrature, charge_n);
  
  ##
  res = zeros (ndofs, 1);
  
  ## ASSEMBLING FIRST ROW
  res(intdofsV)  = A11(intnodes, :) * V;
  res(intdofsV) += A12(intnodes, :) * n;
  
  ## ENFORCING BCs ON FIRST ROW
  for ii = 1:numcontacts
    ddofsV = indexing.V(dnodes{ii});
    
    if (ii == 1) # Metal/semic. interface.
      res(ddofsV) += M(dnodes{ii}, dnodes{ii}) * ...
                     ( (V (dnodes{ii}) - F (device.pins(ii))) - ...
                       (material.PhiB) );
    else # Gate contact.
      res(ddofsV) += M(dnodes{ii}, dnodes{ii}) * ...
                     ( (V (dnodes{ii}) - F (device.pins(ii))) - ... 
                       (material.PhiB + device.Vshift) );
    endif
  endfor
  
  ## COMPUTING SECOND ROW
  A22 = sparse (nnodes, nnodes);

  A22(device.scnodes, device.scnodes) = ...
  bim2a_advection_diffusion ...
    (device.redmsh, mobility.n(! device.insulator) * constants.Vth,
     1, alpha.n(device.scnodes), V(device.scnodes) / constants.Vth);
  
  A22 += bim2a_reaction (device.msh, ! device.insulator, 1 / deltat);

  r22 = A22 * n;
  r22 += bim2a_rhs (device.msh, ! device.insulator, 1 / deltat) .* (-n0); # Avoid cancellation errors.
  
  ## ASSEMBLING SECOND ROW
  res(intdofsn) = r22(intnodes);

  ## ENFORCING BCs ON SECOND ROW
  for ii = 1:numcontacts
    ddofsn = indexing.n(dnodes{ii});
    
    if (ii == 1) # Metal/semic. interface.
        rho = org_gaussian_charge_n (V (dnodes{ii}), material, constants, quadrature);
        
        nimposed = -rho / constants.q;
        
        res(ddofsn) += M(dnodes{ii}, dnodes{ii}) * (n(dnodes{ii}) - nimposed);
    endif
  endfor
  
  ## ADJUST FOR ZERO INSULATOR CHARGE
  if isfield (device, 'scnodes')
    insn = !device.scnodes;
      res(indexing.n(insn)) = M(insn, insn) * n(insn);
  endif
  
  ## ASSEMBLING THIRD ROW
  res(indexing.F) = A * (F - F0) / deltat + ...
                    C + r * I(:);
  
  ## COMPUTING FOURTH ROW
  res(indexing.I) = I;
  
  for ii = 1:numcontacts
    rr = dnodes{device.pins(ii)};
    
    # Displacement current.
    res(indexing.I(ii)) -= ...
      device.section * sum ((A11(rr,:)) * (V - V0)) / deltat;
    res(indexing.I(ii)) -= ...
      device.section * sum ((A12(rr,:)) * (n - n0)) / deltat;
    
    # Electron current.
    res(indexing.I(ii)) -= ...
      - device.section * constants.q * sum (r22(rr));
  endfor
  
  res(indexing.V) /= algorithm.rowscaling(1);
  res(indexing.n) /= algorithm.rowscaling(2);
  res(indexing.F) /= algorithm.rowscaling(3);
  res(indexing.I) /= algorithm.rowscaling(4);
  
  endfunction
