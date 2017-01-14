## Do a "clear all;" before calling the function more than once with different values of "sigman".
function [] = mis_main2d_CF(PhiB = 0.54, sigman = 2.6, mu0n = [], Vshift = [], Csb = [], t_semic = [], voltage = 35)
    
    folder = ['output/CF_' num2str(voltage) '_PhiB_' num2str(PhiB) '_sigman_' num2str(sigman) '_'];
    
    [device, material, constants, ...
     quadrature, charge_n, ...
     algorithm, Vin, nin] = mis_main_params2d(PhiB, sigman, mu0n, Vshift, Csb, t_semic);

     ## 2D simulation
    eps_ins0 = material.eps_ins;
    Csb0 = device.Csb;
    
    CF = csvread ('data_CF_N2200_12-10-27.csv');
    frequency = [CF(2:end, 1); 2 * logspace(5, 11, 91)(2:end)']; # [Hz]
    eps_ins_r = CF(2:end, 5); # [~]

    save ('-mat', [folder 'mis_output_CF_freq.mat'], 'frequency');
    
    loop = tic ();
    for i = 1:numel (frequency)
        freq = frequency(i);
        
        # Compute frequency-dependent eps_ins and Csb.
        if (freq <= frequency(61))
            material.eps_ins_r = interp1(frequency(1:61), eps_ins_r, freq, 'linear', 'extrap'); # [~]
        else
            material.eps_ins_r = real(2.56 + (eps_ins_r(1) - 2.56) ./ (1 + sqrt(-1) * freq * 1e-2).^0.24); # [~]
        end
        material.eps_ins = constants.eps0 * material.eps_ins_r; # [F / m]
        device.Csb = Csb0 * material.eps_ins / eps_ins0; # [F]
        
        # Time span for simulation
        tmin  = -100;
        tmax  = 5 / freq;
        tspan = [tmin -90 -50 linspace(0, tmax, 1000)];
        
        # Circuit boundary conditions and initial condition.
        # Full circuit.
        bc = @(t, F = NA) mis_bcs_AC_circuit(t, voltage, freq, device.Csb,
                                             device.Vshift, F);
        Fin = -device.Vshift * [1 0 device.Csb 0]';
        
        simulation = tic ();
        [Vout, nout, F, Itot, times] = ...
        org_secs2d_newton (device, material, constants, quadrature,
                           charge_n, algorithm,
                           Vin, nin, tspan, bc, Fin);
        
        save ('-mat', sprintf ([folder 'mis_output_CF_%gHz.mat'], freq),
              'device', 'Vout', 'nout', 'F', 'Itot',
              'times', 'tmin', 'tmax', 'freq', 'voltage');
        
        printf ("\n\n1D simulation run time = %g\n", toc(simulation));
    endfor

    printf ("\n\n1D loop run time = %g\n\n", toc(loop));
endfunction
