module fftw_mod
    !integer, parameter :: fftw_forward         = -1, fftw_backward                  = 1
    !integer, parameter :: fftw_real_to_complex = -1, fftw_complex_to_real           = 1
    !integer, parameter :: fftw_estimate        =  0, fftw_measure=1,fftw_use_wisdom = 16
    !integer, parameter :: fftw_in_place        =  8
    integer * 8 :: plan1_rk, plan1_kr, plan2_kr, plan2_rk
endmodule fftw_mod