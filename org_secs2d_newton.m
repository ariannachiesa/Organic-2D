function [V, n, F, I, tout] = org_secs2d_newton ...
                              (device, material, constants,
                               quadrature, charge_n,
                               algorithm, Vin, nin,
                               tspan, va, Fin = [], Iin = [])

  [~, strnm] = fileparts (tmpnam);
  lastsaved = nsaves = 0;
  if (isfield (algorithm, 'plotsOff'));
    plotsOff = algorithm.plotsOff;
  else
    plotsOff = 1;
  endif
  
  if (! isfield (algorithm, 'savedata'));
    algorithm.savedata = 0;
  endif
  
  if (! isfield (algorithm, 'maxdtincr'));
    algorithm.maxdtincr = 2;
  endif

  if (! isfield (algorithm, 'nsteps_check'))
    algorithm.nsteps_check = 3;
  endif

  if (! isfield (algorithm, 'maxit_mnewton'))
    algorithm.maxit_mnewton = 10;
  endif

  if (! isfield (algorithm, 'clamping'))
    algorithm.clamping = algorithm.colscaling / 10;
    warning (['No algorithm.clamping field found.'...
              'Using algorithm.colscaling/10 instead'])
  endif

  newton_solves = modified_newton_solves = 0;
  info.dt = [];
  info.tstep = [];
  info.itnewton = [];
  info.itmodnewton = [];
  info.rejmodnewton = [];
  info.tklast = [];
  info.clamplast = [];
  info.incrlast = [];
  info.incr0last = [];
  info.whyfinished = [];
  rejected = 0;
  
  nnodes = columns(device.msh.p);
  assert (numel (unique (device.pins)) == 2);
  Nelements = columns(device.msh.t);

  if (isfield (algorithm, "dt0"))
    dt = algorithm.dt0;
  else
    dt = (tspan(2) - tspan(1)) / 1000;
  endif

  if (isfield (algorithm, 'dtmax'))
    dt = min(dt, algorithm.dtmax);
  else
    algorithm.dtmax = dt * algorithm.maxdtincr ^ 3;
  endif

  if (! isfield (algorithm, 'dtmin'))
    algorithm.dtmin = 1e-7;
  endif

  if ((! isfield (algorithm, 'restartfile')) || isempty(algorithm.restartfile))
    ## INIT STATE FROM SCRATCH
    tout(tstep = 1) = t = tspan (1);

    [V, n] = deal (Vin, nin);

    [A, B, C, r, F] = va (t);
    Nextvars  = numel (F);
    NI  = numel (device.pins);
    I = zeros (NI, 1);
    if !isempty(Fin)
      F = Fin; 
    end
    if !isempty(Iin)
      I = Iin;
    end
    firstfixtstep = 2;
  else
    ## RESTART LOADING FILE
    disp (' ')
    warning (['Restarting the simulation from file '...
              algorithm.restartfile]);
    disp ('      -  Loading ...')
    restart = load (algorithm.restartfile);
    V = [restart.V2 restart.Vold];
    n = [restart.n2 restart.nold];
    F = [restart.F2 restart.Fold];
    I = [restart.I2 restart.Iold];

    tspan = unique([restart.t restart.told tspan(tspan >= restart.t)]);
    tstep = 2;
    tout = [restart.told restart.t];
    Nextvars  = rows (F);
    NI  = rows (I);
    [A, B, C, r] = va (restart.t, F(:,end));
    nsaves = restart.nsaves;
    firstfixtstep = 3;
    dt0 = dt = restart.dt;
    printf ('      -  Restarting at time %.12g s with dt = %.12g\n\n', tout(end), dt)
  endif

  disp (sprintf ('Connecting contact at x = %.6e to circuit pin %d', 
    device.msh.p(1,1),   device.pins(1)));
  disp (sprintf ('Connecting contact at x = %.6e to circuit pin %d', 
    device.msh.p(1,end), device.pins(2)));

  # node ordering
  indexing.V = 1:2:2*nnodes;
  indexing.n = 2:2:2*nnodes;
  indexing.F = (1:Nextvars) + 2 * nnodes;
  indexing.I = (Nextvars + 2 * nnodes) + (1:NI);

  for n_fix_tstep = firstfixtstep:numel (tspan) ## TIME FIXED SPAN
    t = tspan(n_fix_tstep - 1);
    while (t < tspan(n_fix_tstep)) ## TIME STEP
      printf ('--\n');

      t = tout(++tstep) = min (t + dt, tspan(n_fix_tstep)); 
      dt = t - tout(tstep - 1);

      info.dt = [info.dt, dt];
      info.tstep = [info.tstep, tstep];
      info.whyfinished = [info.whyfinished, NA];

      incr0 = 4 * algorithm.maxnpincr;

      [A, B, C, r] = va (t, F(:, end));
      [V0, n0, F0, I0] = org_secs_state_predict ...
                         (device, material, constants, 
                          algorithm, V, n, F, I, tstep, tout);                          
                          
      [V2, n2, F2, I2] = deal (V0, n0, F0, I0);
      
      in = 0;
      reject = false;

      totmn = rejmn = 0;

      while (!reject && (in < algorithm.maxit)) ## NEWTON STEP
        in += 1;
        incrlast = incr1 = NA;

        [V1, n1, F1, I1] = deal (V2, n2, F2, I2);

        [A, B, C, r] = va (t, F2);

        res = org_secs2d_newton_residual (device, material, constants,
                                          quadrature, charge_n,
                                          algorithm, V2, n2, F2, I2,
                                          V(:, tstep-1), n(:, tstep-1),
                                          F(:, tstep-1), I(:, tstep-1),
                                          dt, A, B, C, r, indexing);

        if (in == 1)
          [~, ~, resall] = compute_residual_norm (res, algorithm, indexing);
          
          algorithm.rowscaling .*= (resall + 1);
          
          res = org_secs2d_newton_residual (device, material, constants,
                                            quadrature, charge_n,
                                            algorithm, V2, n2, F2, I2,
                                            V(:, tstep-1), n(:, tstep-1),
                                            F(:, tstep-1), I(:, tstep-1),
                                            dt, A, B, C, r, indexing);
        endif

        [resnrm(in), whichone, resall] = compute_residual_norm ...
                                           (res, algorithm, indexing); 

        jac = org_secs2d_newton_jacobian (device, material, constants, 
                                          quadrature, charge_n,
                                          algorithm, V2, n2, F2,
                                          dt, A, B, C, r, indexing);

        [jacL, jacU, jacP, jacQ, jacR] = lu (jac);
        delta = - jacQ * (jacU \ (jacL \ (jacP * (jacR \ res))));
        newton_solves +=1;
        
        dV = delta (indexing.V) * algorithm.colscaling(1);
        dn = delta (indexing.n) * algorithm.colscaling(2);
        dF = delta (indexing.F) * algorithm.colscaling(3);
        dI = delta (indexing.I) * algorithm.colscaling(4);
        
        [V2, n2, F2, I2, clamp, tauk] = org_secs_safe_increment ...
                                        (V1, n1, F1, I1,
                                         dV, dn, dF, dI,
                                         algorithm, device, constants);
        
        if ((clamp <= 0) || (tauk <= 0))
          reject = true;
          info.whyfinished(end) = -2;
          break
        endif

        [incr0v, incr0n, incr0F, incr0I] = ...
        compute_variation (V0, n0, F0, I0, 
                           V2, n2, F2, I2, 
                           algorithm, device, material,
                           constants, min (clamp, tauk));

        [incr0, whichone] = ...
        max ([incr0v, incr0n, incr0F, incr0I]);
        
        if (incr0 > algorithm.maxnpincr
            && dt > algorithm.dtmin)
          MAXINCR_MSG (tstep, t, in, whichone, incr0, resall, algorithm)
          reject = true;
          info.whyfinished(end) = -3;
          break;
        endif
        
        [incr1v, incr1n, incr1F, incr1I] = ...
        compute_variation (V1, n1, F1, I1, 
                           V2, n2, F2, I2, 
                           algorithm, device, material,
                           constants, min (clamp, tauk));

        [incr1, whichone] = ...
        max ([incr1v, incr1n, incr1F, incr1I]);
        incnrm(in) = incr1;
        inc_clamp(in) = incr1 / clamp;
        incrlast = incr1;
        
        if (resnrm(in) < algorithm.toll)
          CONV_MSG (tstep, in, 1, t, 'res', 
                    whichone, incr1, resall);
          info.whyfinished(end) = 1;
          break;
        endif

        if (in > algorithm.nsteps_check 
            && incnrm(in) > incnrm(in - algorithm.nsteps_check)
            && resnrm(in) > resnrm(in - algorithm.nsteps_check)
            && dt > algorithm.dtmin)
          DIV_MSG (tstep, t, in, whichone, incnrm, incr1, 
                   resall, algorithm.nsteps_check)
          reject = true;
          info.whyfinished(end) = -4;
          break;
        endif

        printf ('incr (%d) = %.8e ', whichone - 1, incr1);
        printf ('residual = [ %.8e %.8e %.8e %.8e ] \n',
                resall(1), resall(2), resall(3), resall(4));

        if (incr1 < algorithm.toll)
          CONV_MSG (tstep, in, 1, t, 'incr', 
                    whichone, incr1, resall);
          info.whyfinished(end) = 2;
          break;
        endif

        ## MODIFIED NEWTON
        [V1, n1, F1, I1] = deal (V2, n2, F2, I2);
        incrk = NA;
        rejectmnewton = convergedmnewton = false;

        for imn = 1:algorithm.maxit_mnewton
          [Vk, nk, Fk, Ik] = deal (V2, n2, F2, I2);

          [A, B, C, r] = va (t, F2);
          res = org_secs2d_newton_residual (device, material, constants,
                                            quadrature, charge_n,
                                            algorithm, V2, n2, F2, I2,
                                            V(:, tstep-1), n(:, tstep-1),
                                            F(:, tstep-1), I(:, tstep-1),
                                            dt, A, B, C, r, indexing);

          [resnrmk(imn), whichone, resall] = ...
          compute_residual_norm (res, algorithm, indexing); 

          delta = - jacQ * (jacU \ (jacL \ (jacP * (jacR \ res))));
          modified_newton_solves +=1;

          dV = delta (indexing.V) * algorithm.colscaling(1);
          dn = delta (indexing.n) * algorithm.colscaling(2);
          dF = delta (indexing.F) * algorithm.colscaling(3);
          dI = delta (indexing.I) * algorithm.colscaling(4);

          [V2, n2, F2, I2, clamp, tauk] = ...
          org_secs_safe_increment (Vk, nk, Fk, Ik,
                                   dV, dn, dF, dI, 
                                   algorithm, device, constants);
          
          if ((tauk <=0) || (clamp <= 0))
            reject = true;
            info.whyfinished(end) = -5;
            break
          endif

          [incrkv, incrkn, incrkF, incrkI] = ...
          compute_variation (Vk, nk, Fk, Ik, 
                             V2, n2, F2, I2, 
                             algorithm, device, material,
                             constants, min (clamp, tauk));

          [incrk, whichone] = ...
          max ([incrkv, incrkn, incrkF, incrkI]);
          incnrmk(imn) = incrk;
          inck_clamp(imn) = incrk / clamp;
          incrlast = incrk;

          if (resnrmk(imn) < algorithm.toll)
            CONV_MSG (tstep, in, imn, t, 'res', 
                      whichone, incrk, resall);
            info.whyfinished(end) = 3;
            rejectmnewton = false;
            convergedmnewton = true;
            break;
          endif

          if (! plotsOff) 
            figure (1)
            subplot(2, 1, 1)
            semilogy (incnrm(1:in), 'bo-;IN;', 
                      inc_clamp(1:in), 'ms-;INc;', 
                      incnrmk(1:imn), 'gd-;IM;',
                      inck_clamp(1:imn), 'cs-;IMc;' 
                      );
            xlim ([0, (max ([15, in, imn]))]);
            subplot(2, 1, 2)
            semilogy (resnrm(1:in), 'rx-;RN;',
                      resnrmk(1:imn), 'm+-;RM;');
            xlim ([0, (max ([15, in, imn]))]);
            drawnow
          endif

          if (imn > algorithm.nsteps_check
              && (resnrmk(imn) > resnrmk(imn - algorithm.nsteps_check)) 
              && (incnrmk(imn) > incnrmk(imn - algorithm.nsteps_check)))
            DIV_MN_MSG (tstep, t, in, imn, whichone, incnrmk,
                        incrk, resall, algorithm.nsteps_check)
            rejectmnewton = true;
            convergedmnewton = false;
            rejmn++;
            break;
          endif

          printf ('incr (%d) = %.8e ', whichone - 1, incrk);
          printf ('residual = [ %.8e %.8e %.8e %.8e ] * \n',
                  resall(1), resall(2), resall(3), resall(4));

          if (incrk < algorithm.toll)
            CONV_MSG (tstep, in, imn, t, 'incr', 
                      whichone, incrk, resall);
            info.whyfinished(end) = 4;
            rejectmnewton = false;
            convergedmnewton = true;
            break;
          endif

        endfor ## MODIFIED NEWTON

        totmn += imn;

        if (reject)
          break
        endif

        if (rejectmnewton == true)
          [V2, n2, F2, I2] = deal (V1, n1, F1, I1);
        elseif (convergedmnewton == true)
          break;
        endif

        if (exist (fullfile (P_tmpdir, 'stop_secs2d'), 'file'))
          break;
        endif

        if (in >= algorithm.maxit)
          printf ...
            (['\nmaximum number of Newton iterations reached, ', ...
              'try reducing timestep...\n']);
          if (dt > algorithm.dtmin)        
            reject = true;
            info.whyfinished(end) = -6;
          endif
        endif
        
      endwhile ## NEWTON STEP

      info.itnewton = [info.itnewton, in];
      info.itmodnewton = [info.itmodnewton, totmn];
      info.rejmodnewton = [info.rejmodnewton, rejmn];
      info.clamplast = [info.clamplast, clamp];
      info.tklast = [info.tklast, tauk];
      info.incrlast = [info.incrlast, incrlast];
      info.incr0last = [info.incr0last, incr0];
      assert (size(info.whyfinished), size(info.incrlast), 0)

      if (reject)

        ++rejected;
        t = tout (--tstep);
        dt *= algorithm.dtcut;

        printf ('\nreverting to time step %d\n', tstep-1);
        printf ('reducing time step:');
        printf ('\tmodel time %.17g s,', t);
        printf ('\tnew dt %g s\n', dt);

      else
        [V(:, tstep), n(:, tstep), F(:, tstep), I(:, tstep)] = ...
                                            deal (V2, n2, F2, I2);

        dtfact = min (.8 * sqrt (algorithm.maxnpincr / incr0),
                      algorithm.maxdtincr);
        dt = max (min (dtfact * dt, algorithm.dtmax), algorithm.dtmin);

        printf (['\n\tt = %g <= %g; dtfact = %g, ' ...
                 'estimate for next time step size: '...
                 'dt = %g (incr0 = %e)\n'],  t, max(tspan), dtfact, dt, incr0);
      endif
      
      if (algorithm.savedata)
        save ('-binary', '-z', sprintf ('secs2d_%s.octbin.gz', strnm), 
            'F', 'I', 'tout', 'algorithm', 'newton_solves',
            'modified_newton_solves', 'info');
      endif

      if (algorithm.savedata && ((tstep - lastsaved) >= 100))
        lastsaved = tstep;
        Vold = V(:, tstep - 1);
        nold = n(:, tstep - 1);
        Fold = F(:, tstep - 1);
        Iold = I(:, tstep - 1);
        told = tout(tstep - 1);
        save ('-binary', '-z', 
              sprintf ('secs1d_%s_%4.4d.octbin.gz', strnm, ++nsaves), 
              'Vold', 'nold', 'Fold', 'Iold', 'told',
              'V2', 'n2', 'F2', 'I2', 't', 'dt',
              'nsaves', 'algorithm');
        clear Vold nold Fold Iold told
      endif

      if (exist (fullfile (P_tmpdir, 'stop_secs'), 'file'))
        warning ('Stop file detected \n Saving and going off');
        break;
      endif
    endwhile ## TIME STEP

    if (algorithm.savedata)
      lastsaved = tstep;
      Vold = V(:, tstep - 1);
      nold = n(:, tstep - 1);
      Fold = F(:, tstep - 1);
      Iold = I(:, tstep - 1);
      told = tout(tstep - 1);
      save ('-binary', '-z', 
            sprintf ('secs1d_%s_%4.4d.octbin.gz', strnm, ++nsaves), 
            'Vold', 'nold', 'Fold', 'Iold', 'told',
            'V2', 'n2', 'F2', 'I2', 't', 'dt',
            'nsaves', 'algorithm');
      clear Vold nold Fold Iold told
    endif
  endfor ## TIME FIXED SPAN

  printf ('total number of rejected time steps: %d\n', rejected);

endfunction

function [incrV, incrn, incrF, incrI] = compute_variation ...
                                        (Va, na, Fa, Ia,
                                         Vb, nb, Fb, Ib,
                                         algorithm, device,
                                         material, constants, clamping)

  cl = algorithm.clamping;
  sc = device.scnodes;
  ni = material.ni;
  
  incrV = norm (Vb - Va, inf) / (norm (Va, inf) * clamping + cl(1));
  incrn = constants.Vth * norm (log (nb(sc) ./ na(sc)), inf) / ...
          (norm (Va(sc) - constants.Vth * log (na(sc) ./ ni), inf) * clamping + cl(2));
  incrF = norm (Fb - Fa, inf) / (norm (Fa, inf) * clamping + cl(3));
  incrI = norm (Ib - Ia, inf) / (norm (Ia, inf) * clamping + cl(4));

endfunction

function [resnrm, whichone, resall] = compute_residual_norm ...
                                      (res, alg, idx) 

  resall(1) = norm (res(idx.V), inf);
  resall(2) = norm (res(idx.n), inf);
  resall(3) = norm (res(idx.F), inf);
  resall(4) = norm (res(idx.I), inf);

  [resnrm, whichone] =  max (resall);

endfunction

function CONV_MSG(tstep, Nstep, mNstep, t, reason, field, incr, res) 

  printf ('\nat time step %d\n', tstep - 1);
  printf ('fixed point iteration %d', Nstep -1);
  printf (', modified Newton iteration %d', mNstep -1);
  printf (', model time %.17g', t);
  printf (': convergence reached (%s) \n', reason);

  printf ('incr (%d', field - 1);
  printf (') = %.8e residual = [ ' , incr);
  printf ('%.8e ', res(1:end-1));
  printf ('%.8e ]\n\n', res(end));

endfunction

function BADGUESS_MSG (tstep, t, code)
  printf ('\n at time step %d ', tstep - 1)
  printf ('model time %g\n', t);
  printf ('the extrapolation');
  printf (' produced negative concentrations');
  printf ('. Error code %d\n', code);

endfunction

function MAXINCR_MSG (tstep, t, Nstep, field, incr, res, alg)
  
  printf ('\nat time step %d ', tstep - 1);
  printf ('model time %g\n', t);
  printf ('fixed point iteration %d', Nstep -1);
  printf (' the increment in field %d', field - 1);
  printf (' has grown too large: ');
  printf ('%g > %g\n', incr, alg.maxnpincr);

  printf ('incr0 (%d', field - 1);
  printf (') = %.8e residual = [ ' , incr);
  printf ('%.8e ', res(1:end-1));
  printf ('%.8e ]\n', res(end));

endfunction

function DIV_MSG (tstep, t, Nstep, field, incrhist, ...
                  incr, res, nsteps_check)

  printf ('\nat time step %d ', tstep - 1);
  printf ('model time %g\n', t);
  printf ('fixed point iteration %d', Nstep);
  printf (', the Newton algorithm is diverging: ');
  printf ('the increment in field %d', field - 1);
  printf (' is not decreasing : ');
  printf ('%g > %g\n', incrhist(Nstep), 
          incrhist(Nstep - nsteps_check));

  printf ('incr (%d', field - 1);
  printf (') = %.8e residual = [ ' , incr);
  printf ('%.8e ', res(1:end-1));
  printf ('%.8e ]\n', res(end));

endfunction

function DIV_MN_MSG (tstep, t, Nstep, mNstep, field, incrhist, ...
                     incr, res, nsteps_check)

  printf ('\nat time step %d ', tstep - 1);
  printf ('model time %g\n', t);
  printf ('fixed point iteration %d', Nstep);
  printf (' modified Newton iteration %d', mNstep);
  printf (', the Modified Newton algorithm is diverging: ');
  printf ('the increment in field %d', field - 1);
  printf (' is not decreasing : ');
  printf ('%g > %g\n', incrhist(mNstep), 
          incrhist(mNstep - nsteps_check));

  printf ('incr (%d', field - 1);
  printf (') = %.8e residual = [ ' , incr);
  printf ('%.8e ', res(1:end-1));
  printf ('%.8e ]\n', res(end));

endfunction
