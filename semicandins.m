function [scnodes, insel] = ...
         semicandins (imsh, insregions, scregions)
  
  SUBMESHFUN = @msh2m_submesh;

  [~, tmp] = SUBMESHFUN (imsh, [], scregions);
  scnodes = false (columns (imsh.p), 1);
  scnodes (tmp) = true;

  if (! isempty (insregions))
    [~, ~, tmp] = SUBMESHFUN (imsh, [], insregions);
  else
    tmp = [];
  endif
  
  insel = false (columns (imsh.t), 1);
  insel (tmp) = true;

endfunction