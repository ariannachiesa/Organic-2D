##
function [device, material, constants, ...
          quadrature, charge_n, ...
          algorithm, Vin, nin] = mis_main_params2d(PhiB, sigman, mu0n = [], Vshift = [], Csb = [], t_semic = [])
    
    addpath("./bim");
    addpath("./msh");
   
    if (isempty(mu0n))
        mu0n = 4.29110133911508e-6; # [m^2 V^{-1} s^{-1}]
    endif
    
    if (isempty(Vshift))
        Vshift = 1.79738; # [V]
    endif
    
    if (isempty(Csb))
        Csb = 1.16183675549126e-11; # [F]
    endif
    
    if (isempty(t_semic))
        t_semic = 3.49436549222355e-8; # [m]
    endif
    
    ## Physical constants and parameters
    constants = org_secs_physical_constants_fun (295);
    material  = org_secs_N2200_material_properties (constants, PhiB, sigman, mu0n);
    
    device.Vshift = Vshift; # [V]
    device.Csb = Csb; # [F]
    
    t_ins   = 4.41e-7; # [m]
    L       = 1e-5; # [m]
    
    nNodes_x = 50;
    nNodes_y = 20;
    
    nNodesSemic_y = floor (0.6 * nNodes_y);
    nNodesIns_y = nNodes_y - nNodesSemic_y;
    
    x = linspace(0, L, nNodes_x);
        
    ySemic = linspace(-t_semic, 0, nNodesSemic_y);
    yIns   = linspace(0, t_ins, nNodesIns_y + 1);
    
    mshsc = msh2m_structured_mesh(x,ySemic,1,[1,2,5,4]);
    mshins = msh2m_structured_mesh(x,yIns,2,[5,2,3,4]);
    
    mesh = msh2m_join_structured_mesh(mshsc,mshins,5,5);
    device.msh = bim2c_mesh_properties(mesh);    
    
    [device.scnodes, device.insulator] = semicandins(device.msh,2,1);
    
    device.redmsh = mesh2d_reduction(device);
    
    # Contacts
    device.pins = [2 1]; # [Bulk, gate].
    device.contacts = [1 7];
    
    device.alldnodes = [];
    device.dnodes = {};
  
    for cc = device.contacts(:)'
      device.dnodes = {device.dnodes{:}, ...
                (bim2c_unknowns_on_side (device.msh, cc))};
      device.alldnodes = union (device.alldnodes, device.dnodes{end});
    endfor
    
    # Device section [m^2]
    device.section = 0.00000081;

    # Device source-drain electric field [V / m]
    Vdrain = 5; # [V]
    device.Efield = Vdrain / L;

    # Quadrature nodes and weights
    [gx, gw] = ghrule (101);
    quadrature.gx = gx;
    quadrature.gw = gw;

    # Initialize gaussian charge rules
    charge_n = @(phi) org_gaussian_charge_n (phi, material, constants, quadrature);

    # Algorithm parameters
    algorithm.ptoll         = 1e-10;
    algorithm.pmaxit        = 1000;
    algorithm.toll          = 1e-4;
    algorithm.maxit         = 5;
    algorithm.nsteps_check  = 3;
    algorithm.maxit_mnewton = 30;
    algorithm.maxnpincr     = 1e-3;
    algorithm.dt0           = 1e-10;
    algorithm.dtcut         = 0.25;
    algorithm.dtmax         = 1;
    algorithm.dtmin         = 1e-12;
    algorithm.maxdtincr     = 2;
    algorithm.plotsOff      = true;
    algorithm.savedata      = 0;

    ## Scaling and clamping coefficients
    algorithm.colscaling = [1, 1, 1, 1];
    algorithm.rowscaling = [1, 1, 1, 1];
    algorithm.clamping   = [1, 1, 1, 1];

    ## Initial guess
    Vguess = material.PhiB * ones (columns(device.msh.p), 1);

    [Vin, nin, res, niter] = org_secs2d_nlpoisson ...
                            (device, material, constants, algorithm, Vguess, charge_n);
    
    ## Other physical parameters
    material.ni = mean(nin(device.scnodes));
endfunction
