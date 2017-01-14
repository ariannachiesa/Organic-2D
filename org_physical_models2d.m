## Physical models for disordered materials.
function [mobility, alpha, der] = ...
            org_physical_models2d (n, material, device, constants,
                                 quadrature, charge_n)
  vpe = 3;
  nm = vpe ./ sum(1 ./ n(device.msh.t(1:vpe,:)), 1)(:);

  # Mobility model.
  mobility.n = org_secs_mobility_EGDM (material.mu0n, nm,
                                       material.N0, material.sigman_kT,
                                       device, constants, charge_n);
  
  # Generalized Einstein relation.
  alpha.n = org_einstein_n (n, material, constants,
                            quadrature, charge_n); # Piece-wise linear.
  
  
  
  ## Functional derivatives.
  if (nargout == 3)
    # Generalized Einstein relation.
    [~, der.dalpha.n] = org_einstein_n (nm, material, constants,
                                        quadrature, charge_n); # Element-wise constant.
  end
  
endfunction