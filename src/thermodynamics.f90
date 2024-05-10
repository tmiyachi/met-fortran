module thermodynamics
  use iso_fortran_env
  implicit none

  real(real32), parameter :: AIR_R_d = 287. !! 乾燥空気の気体定数 [J/K/Kg]
  real(real32), parameter :: AIR_Cp_d = 1004. !! 乾燥空気の低圧比熱 [J/K/Kg]
  real(real32), parameter :: AIR_Cv_d = 718. !! 乾燥空気の定積比熱 [J/K/Kg]

  real(real32), parameter :: WATER_R_v = 461.50 !! 水蒸気の気体定数 [J/K/Kg]
  real(real32), parameter :: WATER_Cp_v = 1885. !! 水蒸気の低圧比熱 [J/K/Kg]
  real(real32), parameter :: WATER_Cv_v = 1390. !! 水蒸気の定積比熱 [J/K/Kg]

  real(real32), parameter :: WATER_Cp_l = 4186. !! 液体の水の低圧比熱 [J/K/Kg]
  real(real32), parameter :: WATER_Cp_i = 2106. !! 氷の低圧比熱 [J/K/Kg]

  real(real32), parameter :: WATER_Es_0c = 611.2 !! 0度での飽和水蒸気圧 [Pa]

  real(real32), parameter :: WATER_Lv_0c = 2.50084e6 !! 0度での水の蒸発熱 [J/Kg]
  real(real32), parameter :: WATER_Ls_0c = 2.836e6 !! 0度での水の昇華熱 [J/Kg]
  real(real32), parameter :: WATER_Lf_0c = 3.337e5 !! 0度での水の融解熱 [J/Kg]

  real(real32), parameter :: AIR_KAPPA_d = AIR_R_d/AIR_Cp_d !! Rd/Cp of dry air
  real(real32), parameter :: AIR_GAMMA_d = AIR_Cp_d/AIR_Cv_d !! Cp/Cv of dry air
  real(real32), parameter :: AIR_EPSILON_d = AIR_R_d/WATER_R_v !! Rd/Rv

  private
  public :: equivalent_potential_temperature, &
            mixing_ratio, &
            mixing_ratio_from_relative_humidity, &
            mixing_ratio_from_specific_humidity, &
            mixing_ratio_from_dewpoint, &
            moist_adiabatic_lapse_rate, &
            lifting_condensation_level, &
            potential_temerature, &
            relative_humidity_from_dewpoint, &
            relative_humidity_from_mixing_ratio, &
            relative_humidity_from_specific_humidty, &
            saturation_equivalent_potential_temperature, &
            saturation_mixing_ratio, &
            saturation_vapor_pressure, &
            specific_humidity_from_dewpoint, &
            specific_humidity_from_mixing_ratio, &
            virtual_potential_temperature, &
            virtual_temperature, &
            virtual_temperature_from_dewpoint, &
            wet_bulb_potential_temperature_davies_jones, &
            wet_bulb_temperature_stull

contains

  function equivalent_potential_temperature(p, t, td) result(theta_e)
    !! Bolton(1980) の近似式を用いて相当温位を計算
    !!
    !! ### Formulation
    !! Bolton(1980) の (15),(24),(39) 式
    !!
    !! $$\theta_{DL} = T\left(\frac{1000}{p - e}\right) ^{R_d/C_pd}$$
    !! $$T_L = \frac{1}{\frac{1}{T_d - 56}+\frac{\ln(T/T_d) }{800}}$$
    !! $$\theta_e = \theta_{DL}\exp\left(\frac{3036}{T_L}-1.78\right) r(1 + 0.448r) $$
    !!
    !! ### Reference
    !! Bolton, D., 1980:The computation of equivalent potential temperature.Mon.Wea.Rev., 108, 1046 - 1053

    real(real32), intent(in) :: p
      !! 気圧 [hPa]
    real(real32), intent(in) :: t
      !! 気温 [K]
    real(real32), intent(in) :: td
      !! 露点温度 [K]
    real(real32) :: theta_e
      !! 相当温位 [K]
    real(real32) :: w, e, t_l, theta_d, theta_dl

    w = mixing_ratio(td, p) ! [g/g]
    e = saturation_vapor_pressure(td) ! [hPa]

    t_l = 1./(1/(td - 56.) + log(t/td)/800.) + 56. ! B80 Eq.15
    theta_d = potential_temerature(p - e, t)
    theta_dl = theta_d*((t/t_l)**(0.28*w)) ! B80 Eq.24

    theta_e = theta_dl*exp((3036./t_l - 1.78)*w*(1 + 0.448*w)) ! B80 Eq.39

  end function equivalent_potential_temperature

  function mixing_ratio(e, p) result(w)
    !! 混合比を計算
    !!
    !! ### Formulation
    !! $$w=\epsilon\frac{e}{p-e}$$

    real(real32), intent(in) :: e
      !! 水蒸気圧 [hPa]
    real(real32), intent(in) :: p
      !! 気圧 (全圧) [hPa]
    real(real32) :: w
      !! 混合比 [g/g]

    w = AIR_EPSILON_d*e/(p - e) ! BA23 eq.5.14

  end function mixing_ratio

  function mixing_ratio_from_relative_humidity(p, t, rh) result(w)
    !! 混合比を相対湿度から計算
    !!
    !! ### Formulation
    !! WMO(2020) の (4.A.16) 式
    !!
    !! $$RH = \frac{w}{w_s}\frac{\epsilon + w_s}{\epsilon + w} \times 100$$
    !!
    !! ### Reference
    !!   WMO, 2020: Guide to Meteorological Instruments and Methods of Observation, Volume 1:
    !!     Measurement of Meteorological Variables., https://library.wmo.int/idurl/4/68695

    real(real32), intent(in) :: p
      !! 気圧 [hPa]
    real(real32), intent(in) :: t
      !! 気温 [K]
    real(real32), intent(in) :: rh
      !! 相対湿度 [%]
    real(real32) :: w
      !! 混合比 [g/g]
    real(real32) :: r, w_s

    r = rh/100. ! unit less ratio
    w_s = saturation_mixing_ratio(p, t) ! [g/g]
    w = AIR_EPSILON_d*w_s*r/(AIR_EPSILON_d + w_s*(1 - r))

  end function mixing_ratio_from_relative_humidity

  function mixing_ratio_from_specific_humidity(q) result(w)
    !! 混合比を比湿から計算
    !!
    !! ### Formulation
    !! Bohren & Albrecht (2023) の (5.10) 式
    !!
    !! $$w=\frac{q}{1-q}$$
    !!
    !! ### Reference
    !!   Bohren, C. F., and B. A. Albrecht, 2023: Atmospheric Thermodynamics Second Edition.
    !!     Oxford University Press, 579 pp

    real(real32), intent(in) :: q
      !! 比湿 [g/g]
    real(real32) :: w
      !! 混合比 [g/g]

    w = q/(1 - q) ! BA23 eq.5.10

  end function mixing_ratio_from_specific_humidity

  function mixing_ratio_from_dewpoint(p, td) result(w)
    !! 混合比を露点温度から計算

    real(real32), intent(in) :: p
      !! 気圧 [hPa]
    real(real32), intent(in) :: td
      !! 露点温度 [K]
    real(real32) :: w
      !! 混合比 [g/g]
    real(real32) :: e

    e = saturation_vapor_pressure(td) ! [hPa]
    w = mixing_ratio(e, p)

  end function mixing_ratio_from_dewpoint

  function moist_adiabatic_lapse_rate(p, t) result(dtdp)
    !! 飽和空気塊の湿潤断熱減率を Bohren & Albrecht の近似式を用いて計算する
    !!
    !! ### Formulation
    !! Bohren & Albrecht の近似式 (6.116)
    !!
    !! $$ \frac{dT}{dz} = -\frac{g}{c_{pd}}\frac{1 + \frac{l_vw_s}{R_dT}}{1 + \frac{l_v^2w_s}{c_{pd}R_vT^2}}$$
    !!
    !! を静水圧平衡の式で変形する
    !!
    !! ### Reference
    !!   Bohren, C. F., and B. A. Albrecht, 2023: Atmospheric Thermodynamics Second Edition.
    !!     Oxford University Press, 579 pp

    real(real32), intent(in) :: p
      !! 気圧 [hPa]
    real(real32), intent(in) :: t
      !! 気温 [K]
    real(real32) :: dtdp
      !! 湿潤断熱減率 [K/hPa]
    real(real32) :: w_s

    w_s = mixing_ratio_from_dewpoint(p, t) ! [g/g] 飽和空気塊なのでt=td
    dtdp = 1./p*(AIR_R_d*t + WATER_Lv_0c*w_s)/(AIR_Cp_d + (WATER_Lv_0c**2*w_s/WATER_R_v/t/t)) ! BA23 eq.6.116

  end function moist_adiabatic_lapse_rate

  subroutine lifting_condensation_level(p, t, td, p_lcl, t_lcl)

    !! 持ち上げ凝結高度を  Bohren & Albrecht (2003) の近似式を用いて計算
    !!
    !! ### Formulation
    !! Bohren & Albrecht (2023) の (6.32) 式
    !!
    !! $$T_{LCL}=\frac{1-AT_d}{1/T_d + B\ln(T/T_d)-A}$$
    !! $$A=-\left(\frac{c_{pv}-c_{pw}}{l_r}-\frac{c_{pd}}{\epsilon l_v}\right),\,\,B=\frac{c_{pd}}{\epsilon l_v}$$
    !! $$l_v=l_{vr} - (c_{pv}-c_{pw})T_r$$
    !!
    !! ### Reference
    !!   Bohren, C. F., and B. A. Albrecht, 2023: Atmospheric Thermodynamics Second Edition.
    !!     Oxford University Press, 579 pp

    real(real32), intent(in) :: p
    real(real32), intent(in) :: t
    real(real32), intent(in) :: td
    real(real32), intent(out) :: p_lcl
    real(real32), intent(out) :: t_lcl
    ! real(real32), parameter :: a = 1.26e-3, b = 5.14e-4
    real(real32), parameter :: lr = WATER_Lv_0c - (WATER_Cp_v - WATER_Cp_l)*273.15, &
                               a = -(WATER_Cp_v - WATER_Cp_l)/lr*AIR_Cp_d/AIR_KAPPA_d/lr, &
                               b = AIR_Cp_d/AIR_KAPPA_d/lr

    t_lcl = (1 - a*td)/(1/td + b*log(t/td) - a)
    p_lcl = p*(t_lcl/t)**(1./AIR_KAPPA_d)

  end subroutine lifting_condensation_level

  function potential_temerature(p, t) result(theta)
    !! 温位を計算
    !!
    !! ### Forumulation
    !! $$\theta=T\left(\frac{1000}{p}\right)^{R_{d}/c_{pd}}$$

    real(real32) :: p
      !! 気圧 [hPa]
    real(real32) :: t
      !! 気温 [K]
    real(real32) :: theta
      !! 温位 [K]

    theta = t*(1000./p)**AIR_KAPPA_d

  end function potential_temerature

  function relative_humidity_from_dewpoint(t, td) result(rh)
    !! 相対湿度を露点温度から計算
    !!
    !! ### Forumulation
    !! WMO(2020) の定義式 (4.A.15)
    !!
    !! $$RH=\frac{e}{e_s} \times 100$$
    !!
    !! ### Reference
    !!   WMO, 2020: Guide to Meteorological Instruments and Methods of Observation, Volume 1:
    !!     Measurement of Meteorological Variables., https://library.wmo.int/idurl/4/68695

    real(real32) :: t
      !! 気温 [K]
    real(real32) :: td
      !! 露点温度 [K]
    real(real32) :: rh
      !! 相対湿度 [%]

    rh = saturation_vapor_pressure(td)/saturation_vapor_pressure(t)*100.

  end function relative_humidity_from_dewpoint

  function relative_humidity_from_mixing_ratio(p, t, w) result(rh)
    !! 相対湿度を混合比から計算
    !!
    !! ### Forumulation
    !! WMO(2020) の (4.A.16) 式
    !!
    !! $$RH = \frac{w}{w_s}\frac{\epsilon + w_s}{\epsilon + w} \times 100$$
    !!
    !! ### Reference
    !!   WMO, 2020: Guide to Meteorological Instruments and Methods of Observation, Volume 1:
    !!     Measurement of Meteorological Variables., https://library.wmo.int/idurl/4/68695

    real(real32), intent(in) :: p
      !! 気圧 [hPa]
    real(real32), intent(in) :: t
      !! 気温 [K]
    real(real32), intent(in) :: w
      !! 混合比 [g/g]
    real(real32) :: rh
      !! 相対湿度 [%]
    real(real32) :: w_s

    w_s = saturation_mixing_ratio(p, t) ! [g/g]
    rh = w/w_s*(AIR_EPSILON_d + w)*(AIR_EPSILON_d + w)*100

  end function relative_humidity_from_mixing_ratio

  function relative_humidity_from_specific_humidty(p, t, q) result(rh)
    !! 相対湿度を比湿から計算

    real(real32), intent(in) :: p
      !! 気圧 [hPa]
    real(real32), intent(in) :: t
      !! 気温 [K]
    real(real32), intent(in) :: q
      !! 比湿 [g/g]
    real(real32) :: rh
      !! 相対湿度 [%]
    real(real32) :: w

    w = mixing_ratio_from_specific_humidity(q)
    rh = relative_humidity_from_mixing_ratio(p, t, w)

  end function relative_humidity_from_specific_humidty

  function saturation_equivalent_potential_temperature(p, t) result(theta_es)
    !! 飽和相当温位を Bolton(1980) の近似式を用いて計算
    !!
    !! ### Formulation
    !!
    !! Bolton(1980) の (15),(24),(39) 式で空気塊が飽和している \(T=T_l\), \(e=e_s\) として計算
    !!
    !! $$\theta_{DL} = T\left(\frac{1000}{p - e_s}\right) ^{R_d/C_pd}$$
    !! $$\theta_e = \theta_{DL}\exp\left(\frac{3036}{T}-1.78\right) r(1 + 0.448r) $$
    !!
    !! ### Reference:
    !!   Bolton, D., 1980: The computation of equivalent potential temperature.
    !!     Mon. Wea. Rev., 108, 1046-1053

    real(real32), intent(in) :: p
      !! 気圧 [hPa]
    real(real32), intent(in) :: t
      !! 気温 [K]
    real(real32) :: theta_es
      !! 飽和相当温位 [K]
    real(real32) :: w, e, t_l, theta_dl

    w = saturation_mixing_ratio(p, t) ! [g/g]
    e = saturation_vapor_pressure(t) ! [hPa]

    theta_dl = potential_temerature(p - e, t) ! B80 Eq.24 when T_k=T_l
    theta_es = theta_dl*exp((3036./t_l - 1.78)*w*(1 + 0.448*w)) ! B80 Eq.39 when T_l = T_k

  end function saturation_equivalent_potential_temperature

  function saturation_mixing_ratio(p, t) result(w_s)
    !! Bolton(1980) の近似式から飽和混合比を計算

    real(real32), intent(in) :: p
      !! 気圧 [hPa]
    real(real32), intent(in) :: t
      !! 気温 [K]
    real(real32) :: w_s
      !! 飽和混合比 [g/g]
    real(real32) :: e_s

    e_s = saturation_vapor_pressure(t)
    w_s = mixing_ratio(e_s, p)

  end function saturation_mixing_ratio

  function saturation_vapor_pressure(t) result(e_s)
    !! Bolton(1980) の近似式から水に対する飽和水蒸気圧を計算
    !!
    !! ### Formulation
    !! Bolton(1980) の (10) 式
    !!
    !! $$e_s(T) = e_{s0} \exp\left(\frac{17.67T}{T+243.5}\right)$$
    !!
    !! ### Reference
    !! Bolton, D., 1980:The computation of equivalent potential temperature.Mon.Wea.Rev., 108, 1046 - 1053

    real(real32), intent(in) :: t !! Temperature [K]
    real(real32) :: e_s !! Satiration Vapor Pressure [hPa]
    real(real32) :: t_c

    t_c = t - 273.15 ! K -> C
    e_s = WATER_Es_0c*exp(17.67*t_c/(t_c + 243.5)) ! B80 eq.10

  end function saturation_vapor_pressure

  function specific_humidity_from_dewpoint(p, td) result(q)
    !! 比湿を露点温度から計算

    real(real32), intent(in) :: p
      !! 気圧 [hPa]
    real(real32), intent(in) :: td
      !! 露点温度 [K]
    real(real32) :: q
      !! 比湿 [g/g]
    real(real32) :: w

    w = mixing_ratio_from_dewpoint(p, td)
    q = specific_humidity_from_mixing_ratio(w)

  end function specific_humidity_from_dewpoint

  function specific_humidity_from_mixing_ratio(w) result(q)
    !! 比湿を混合比から計算
    !!
    !! ### Formulation
    !! Bohren & Albrecht (2023) の (5.10) 式
    !!
    !! $$q=\frac{w}{1+w}$$
    !!
    !! ### Reference
    !!   Bohren, C. F., and B. A. Albrecht, 2023: Atmospheric Thermodynamics Second Edition.
    !!     Oxford University Press, 579 pp

    real(real32), intent(in) :: w
      !! 混合比 [g/g]
    real(real32) :: q
      !! 比湿 [g/g]

    q = w/(1 + w)

  end function specific_humidity_from_mixing_ratio

  function virtual_potential_temperature(p, t, w) result(theta_v)
    !! 仮温位を計算
    !!
    !! ###Formulation
    !! $$\theta_v=\theta\frac{w+\epsilon}{\epsilon(1+w)}$$

    real(real32), intent(in) :: p
      !! 気圧 [hPa]
    real(real32), intent(in) :: t
      !! 気温 [K]
    real(real32), intent(in) :: w
      !! 混合比 [g/g]
    real(real32) :: theta_v
      !! 仮温位 [K]
    real(real32) :: theta

    theta = potential_temerature(p, t)
    theta_v = theta*(w + AIR_EPSILON_d)/AIR_EPSILON_d/(1 + w)

  end function virtual_potential_temperature

  function virtual_temperature(t, w) result(t_v)
    !! 仮温度を計算
    !!
    !! ###Formulation
    !! $$T_v=T\frac{w+\epsilon}{\epsilon(1+w)}$$

    real(real32), intent(in) :: t
      !! 気温 [K]
    real(real32), intent(in) :: w
      !! 混合比 [g/g]
    real(real32) :: t_v
      !! 仮温度 [K]

    t_v = t*(w + AIR_EPSILON_d)/AIR_EPSILON_d/(1 + w)

  end function virtual_temperature

  function virtual_temperature_from_dewpoint(p, t, td) result(t_v)
    !! 仮温度を露点温度を用いて計算

    real(real32), intent(in) :: p
      !! 気圧 [hPa]
    real(real32), intent(in) :: t
      !! 気温 [K]
    real(real32), intent(in) :: td
      !! 露点温度
    real(real32) :: t_v
      !! 仮温度 [K]
    real(real32) :: w

    w = mixing_ratio_from_dewpoint(p, td)
    t_v = virtual_temperature(t, w)

  end function virtual_temperature_from_dewpoint

  function wet_bulb_potential_temperature_davies_jones(p, t, td) result(theta_w)
    !! 湿球温位を Davies-Jones(2008) の近似式を用いて計算
    !!
    !! ###Formulation
    !! Davies-Jones(2008) の (3.8) 式
    !!
    !! \begin{equation*}
    !! \theta_w =
    !!   \begin{cases}
    !!     \theta_e - \exp\left(\frac{a_0+a_1X+a_2X^2+a_3X^3+a_4X^4}{1+b_1X+b_2X^2+b_3X^3+b_4X^4}\right) & \theta_e \ge 173.15\text{K} \\
    !!     \theta_e & \theta_e < 173.15\text{K}
    !!   \end{cases}
    !! \end{equation*}
    !!
    !! $$X=\frac{\theta_e}{273.15}$$
    !!
    !! ###Reference
    !!   Davies-Jones, R., 2008: An Efficient and Accurate Method for Computing the Wet-Bulb Temperature
    !!     along Pseudoadiabats. Mon. Wea. Rev., 136, 2764-2785, https://doi.org/10.1175/2007MWR2224.1

    real(real32), intent(in) :: p
      !! 気圧 [hPa]
    real(real32), intent(in) :: t
      !! 気温 [K]
    real(real32), intent(in) :: td
      !! 露点温度 [K]
    real(real32) :: theta_w
      !! 湿球温位 [K]
    real(real32) :: theta_e, x, x2, x3, x4, a, b

    theta_e = equivalent_potential_temperature(p, t, td)
    x = theta_e/273.15
    x2 = x*x
    x3 = x2*x
    x4 = x3*x

    if (theta_e >= 173.15) then
      a = 7.101574 - 20.68208*x + 16.11182*x2 + 2.574631*x3 - 5.205688*x4
      b = 1 - 3.552497*x + 3.781782*x2 - 0.6899655*x3 - 0.5929340*x4
      theta_w = theta_e - exp(a/b) ! DJ08 eq.3.8
    else
      theta_w = theta_e
    end if

  end function wet_bulb_potential_temperature_davies_jones

  function wet_bulb_temperature_stull(t, rh) result(tw)
    !! 湿球温度を Stull(2011) の経験式を用いて計算
    !!
    !! 適用可能範囲は原著論文を確認して十分注意すること.
    !!
    !! ###Formulation
    !!
    !! Stull(2011) の (2) 式
    !!
    !! \begin{equation*}
    !! T_w=T\arctan(0.151977(RH+8.313659)^{1/2}+\arctan(T+RH)-\arctan(RH-1.676331) \\
    !!    +0.00391838*(RH)^{3/2}*\arctan(0.023101*RH)-4.686035
    !! \end{equation*}
    !!
    !! ###Refference
    !!   Stull, R., 2011: Wet-Bulb Temperature from Relative Humidity and Air Temperature.
    !!     J.Appl.Meteor.Climatol., 50, 2267 - 2269, https://doi.org/10.1175/JAMC-D-11-0143.1

    real(real32), intent(in) :: t
      !! 気温 [K]
    real(real32), intent(in) :: rh
      !! 相対湿度 [%]
    real(real32) :: tw
      !! 湿球温度 [K]

    tw = (t - 273.15)*atan(0.151977*sqrt(rh + 8.313659)) + atan(t - 273.15 + rh) - atan(rh - 1.676331) &
         + 0.00391838*(rh**1.5)*atan(0.023101*rh) - 4.686035 + 273.15 ! S11 eq.1

  end function wet_bulb_temperature_stull

end module thermodynamics
