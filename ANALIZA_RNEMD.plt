reset
set encoding utf8

material   = "Ti3C2"
potential  = "IFF"
anisotropy = "0"
Lx_target_nm = "50"
Ly_target_nm = "5"
tag = "Ti3C2_IFF_Hu0deg_2layers50x5nm1test"
logfile = "log.lammps"
profile = "profilTi3C2_IFF_Hu0deg_2layers50x5nm1test.txt"

nvt_start   = 0
nvt_end     = 200000
npt_start   = 200000
npt_end     = 400000
rnemd_start = 400000
rnemd_end   = 2400000 

h_ps = 0.001

# f_heat -> eV; units real: kcal/mol 0.0433641153 , units metal: 1.0
heat_to_eV = 0.0433641153

avg_start = 1300000
avg_end   = 2400000

T_target = 300.0
P_target = 0.0
d = 14.6300320000

Nbin_rnemd = 100
swap_every = 200

# 1: NPT = Step PE KE Etot Temp Press Pxx Pyy Lx Ly
# 0: stary log bez Pxx/Pyy/Lx/Ly
npt_has_dimensions = 1

# opcjonalnie: Etot podczas RNEMD
# jesli thermo RNEMD zawiera np. Step Temp Press f_heat Xlo Xhi Ylo Yhi Etot, 1 i podanie numeru kolumny
rnemd_has_etot = 1
rnemd_etot_col = 10

# statystyka
block_ps      = 100.0
conv_every_ps = 10.0
stab_fraction = 0.50

# fit profilu T; 0 = zakresy wpisane recznie, 1 = automatyczne
# dla 20 plastrow auto daje L: 3-9, R: 13-19
auto_fit_range = 1
fit_L_first_bin = 3
fit_L_last_bin  = 9
fit_R_first_bin = 13
fit_R_last_bin  = 19

# kryteria robocze QC
T_mean_tol_percent        = 2.0
T_drift_tol_percent       = 1.0
PE_drift_sigma_max        = 0.50
Etot_drift_sigma_max      = 0.50
geom_drift_tol_percent    = 0.10
heat_R2_min               = 0.999
profile_R2_min            = 0.98
gradient_asym_max_percent = 10.0
kappa_SEM_max_percent     = 10.0
kappa_trend_max_percent   = 5.0
kappa_block_agreement_max_percent = 5.0
deltaT_rel_warn_percent   = 20.0
min_blocks_for_stats      = 3
min_profiles_per_block    = 5

# format logu
# NVT:   Step PE KE Etot Temp Press
# NPT:   Step PE KE Etot Temp Press Pxx Pyy Lx Ly
# RNEMD: Step Temp Press f_heat Xlo Xhi Ylo Yhi [Etot]

# CZASY I KONTROLA WEJSCIA

Nsteps_nvt   = nvt_end-nvt_start
Nsteps_npt   = npt_end-npt_start
Nsteps_rnemd = rnemd_end-rnemd_start
Nsteps_avg   = avg_end-avg_start

t_nvt_ps   = Nsteps_nvt*h_ps
t_npt_ps   = Nsteps_npt*h_ps
t_rnemd_ps = Nsteps_rnemd*h_ps
t_avg_ps   = Nsteps_avg*h_ps

block_steps = int(block_ps/h_ps + 0.5)
conv_every  = int(conv_every_ps/h_ps + 0.5)

if (avg_start < rnemd_start || avg_end > rnemd_end || avg_end <= avg_start) {
    print "BLAD: avg_start/avg_end musi lezec wewnatrz RNEMD."
    exit
}
if (block_steps <= 0 || conv_every <= 0) {
    print "BLAD: block_ps i conv_every_ps musza byc > 0."
    exit
}
if (Nbin_rnemd <= 0 || int(Nbin_rnemd/2)*2 != Nbin_rnemd) {
    print "BLAD: Nbin_rnemd musi byc dodatnia liczba parzysta."
    exit
}
if (block_steps > Nsteps_avg) {
    print "UWAGA: block_ps jest dluzszy niz zakres steady-state; SEM blokowe nie bedzie dostepne."
}
if (conv_every >= Nsteps_avg) {
    conv_every = int(Nsteps_avg/4.0)
    if (conv_every < 1) {conv_every = 1}
    conv_every_ps = conv_every*h_ps
    print "UWAGA: conv_every_ps zmniejszono automatycznie do ",conv_every_ps," ps."
}

nvt_stab_start = nvt_end - stab_fraction*Nsteps_nvt
npt_stab_start = npt_end - stab_fraction*Nsteps_npt
nvt_stab_time_ps = (nvt_end-nvt_stab_start)*h_ps
npt_stab_time_ps = (npt_end-npt_stab_start)*h_ps

# STABILNOSC NVT

stats logfile using \
((valid(6) && $3>0 && $1>=nvt_stab_start && $1<nvt_end) ? ($1-nvt_start)*h_ps : 1/0):\
((valid(6) && $3>0 && $1>=nvt_stab_start && $1<nvt_end) ? $5 : 1/0) \
name "NVT_T" nooutput

stats logfile using \
((valid(6) && $3>0 && $1>=nvt_stab_start && $1<nvt_end) ? ($1-nvt_start)*h_ps : 1/0):\
((valid(6) && $3>0 && $1>=nvt_stab_start && $1<nvt_end) ? $2 : 1/0) \
name "NVT_PE" nooutput

nvt_T_mean  = NVT_T_mean_y
nvt_T_std   = NVT_T_stddev_y
nvt_T_slope = NVT_T_slope
nvt_PE_mean  = NVT_PE_mean_y
nvt_PE_std   = NVT_PE_stddev_y
nvt_PE_slope = NVT_PE_slope

nvt_T_dev_percent   = 100.0*abs(nvt_T_mean-T_target)/T_target
nvt_T_drift_percent = 100.0*abs(nvt_T_slope)*nvt_stab_time_ps/T_target
nvt_PE_shift_sigma  = (nvt_PE_std > 0) ? abs(nvt_PE_slope)*nvt_stab_time_ps/nvt_PE_std : 999.0

nvt_ok = (nvt_T_dev_percent <= T_mean_tol_percent && \
          nvt_T_drift_percent <= T_drift_tol_percent && \
          nvt_PE_shift_sigma <= PE_drift_sigma_max)

# STABILNOSC NPT

stats logfile using \
((valid(6) && $3>0 && $1>npt_stab_start && $1<npt_end) ? ($1-npt_start)*h_ps : 1/0):\
((valid(6) && $3>0 && $1>npt_stab_start && $1<npt_end) ? $5 : 1/0) \
name "NPT_T" nooutput

stats logfile using \
((valid(6) && $3>0 && $1>npt_stab_start && $1<npt_end) ? $6 : 1/0) \
name "NPT_P" nooutput

npt_T_mean  = NPT_T_mean_y
npt_T_std   = NPT_T_stddev_y
npt_T_slope = NPT_T_slope
npt_P_mean  = NPT_P_mean
npt_P_std   = NPT_P_stddev

npt_T_dev_percent   = 100.0*abs(npt_T_mean-T_target)/T_target
npt_T_drift_percent = 100.0*abs(npt_T_slope)*npt_stab_time_ps/T_target
npt_thermal_ok = (npt_T_dev_percent <= T_mean_tol_percent && \
                  npt_T_drift_percent <= T_drift_tol_percent)

npt_Pxx_mean = 0.0
npt_Pyy_mean = 0.0
npt_Pxx_std  = 0.0
npt_Pyy_std  = 0.0
npt_Lx_mean  = 0.0
npt_Ly_mean  = 0.0
npt_Lx_std   = 0.0
npt_Ly_std   = 0.0
npt_Lx_drift_percent = 0.0
npt_Ly_drift_percent = 0.0
npt_geom_ok = 1

if (npt_has_dimensions == 1) {
    stats logfile using \
    ((valid(10) && $3>0 && $1>npt_stab_start && $1<npt_end) ? $7 : 1/0) \
    name "NPT_PXX" nooutput

    stats logfile using \
    ((valid(10) && $3>0 && $1>npt_stab_start && $1<npt_end) ? $8 : 1/0) \
    name "NPT_PYY" nooutput

    stats logfile using \
    ((valid(10) && $3>0 && $1>npt_stab_start && $1<npt_end) ? ($1-npt_start)*h_ps : 1/0):\
    ((valid(10) && $3>0 && $1>npt_stab_start && $1<npt_end) ? $9 : 1/0) \
    name "NPT_LX" nooutput

    stats logfile using \
    ((valid(10) && $3>0 && $1>npt_stab_start && $1<npt_end) ? ($1-npt_start)*h_ps : 1/0):\
    ((valid(10) && $3>0 && $1>npt_stab_start && $1<npt_end) ? $10 : 1/0) \
    name "NPT_LY" nooutput

    npt_Pxx_mean = NPT_PXX_mean
    npt_Pyy_mean = NPT_PYY_mean
    npt_Pxx_std  = NPT_PXX_stddev
    npt_Pyy_std  = NPT_PYY_stddev

    npt_Lx_mean = NPT_LX_mean_y
    npt_Ly_mean = NPT_LY_mean_y
    npt_Lx_std  = NPT_LX_stddev_y
    npt_Ly_std  = NPT_LY_stddev_y

    npt_Lx_drift_percent = 100.0*abs(NPT_LX_slope)*npt_stab_time_ps/npt_Lx_mean
    npt_Ly_drift_percent = 100.0*abs(NPT_LY_slope)*npt_stab_time_ps/npt_Ly_mean

    npt_geom_ok = (npt_Lx_drift_percent <= geom_drift_tol_percent && \
                   npt_Ly_drift_percent <= geom_drift_tol_percent)
}

npt_ok = (npt_thermal_ok && npt_geom_ok)

# RNEMD - TEMPERATURA I GEOMETRIA

# > rnemd_start usuwa wspolny wiersz graniczny NPT/RNEMD
stats logfile using \
((valid(8) && $1>rnemd_start && $1<=rnemd_end) ? ($1-rnemd_start)*h_ps : 1/0):\
((valid(8) && $1>rnemd_start && $1<=rnemd_end) ? $3 : 1/0) \
name "TRN" nooutput

stats logfile using \
((valid(8) && $1>=avg_start && $1<=avg_end) ? ($1-rnemd_start)*h_ps : 1/0):\
((valid(8) && $1>=avg_start && $1<=avg_end) ? $3 : 1/0) \
name "TAVG" nooutput

stats logfile using \
((valid(8) && $1>rnemd_start && $1<=rnemd_end) ? ($7-$6) : 1/0) \
name "LXRN" nooutput

stats logfile using \
((valid(8) && $1>rnemd_start && $1<=rnemd_end) ? ($9-$8) : 1/0) \
name "LYRN" nooutput

stats logfile using ((valid(8) && $1>=avg_start && $1<=avg_end) ? $6 : 1/0) name "XLO" nooutput
stats logfile using ((valid(8) && $1>=avg_start && $1<=avg_end) ? $7 : 1/0) name "XHI" nooutput
stats logfile using ((valid(8) && $1>=avg_start && $1<=avg_end) ? $8 : 1/0) name "YLO" nooutput
stats logfile using ((valid(8) && $1>=avg_start && $1<=avg_end) ? $9 : 1/0) name "YHI" nooutput

Lx = XHI_mean-XLO_mean
Ly = YHI_mean-YLO_mean
A_A2 = Ly*d
Lx_nm = Lx/10.0
Ly_nm = Ly/10.0

rnemd_T_mean  = TAVG_mean_y
rnemd_T_std   = TAVG_stddev_y
rnemd_T_slope = TAVG_slope
rnemd_T_dev_percent   = 100.0*abs(rnemd_T_mean-T_target)/T_target
rnemd_T_drift_percent = 100.0*abs(rnemd_T_slope)*t_avg_ps/T_target
rnemd_temp_ok = (rnemd_T_dev_percent <= T_mean_tol_percent && \
                 rnemd_T_drift_percent <= T_drift_tol_percent)

rnemd_E_shift_sigma = 0.0
energy_ok = 1
if (rnemd_has_etot == 1) {
    stats logfile using \
    (($1>=avg_start && $1<=avg_end) ? ($1-avg_start)*h_ps : 1/0):\
    (($1>=avg_start && $1<=avg_end) ? column(rnemd_etot_col) : 1/0) \
    name "RETOT" nooutput

    rnemd_E_shift_sigma = (RETOT_stddev_y > 0) ? abs(RETOT_slope)*t_avg_ps/RETOT_stddev_y : 999.0
    energy_ok = (rnemd_E_shift_sigma <= Etot_drift_sigma_max)
}

# PROFIL TEMPERATURY

# profil fix ave/chunk: naglowek numeryczny ma $2=Nchunks>1,
# a wiersze chunkow maja Coord1 w jednostkach reduced: 0..1
stats profile using (($2>1.0) ? $1 : 1/0) name "PSTEP" nooutput
stats profile using (($2>1.0) ? $2 : 1/0) name "CHUNKS" nooutput
stats profile using (($2>1.0 && $1>=avg_start && $1<=avg_end) ? $1 : 1/0) name "PAVG" nooutput

profile_first = PSTEP_min
profile_last  = PSTEP_max
nprofiles     = PSTEP_records
nprofiles_avg = PAVG_records
nchunks       = int(CHUNKS_mean)

if (nprofiles < 2) {
    print "BLAD: profil temperatury ma mniej niz dwa zapisy."
    exit
}

profile_every = (profile_last-profile_first)/(nprofiles-1)
dx = Lx/nchunks
profiles_per_block = block_steps/profile_every

if (nchunks != Nbin_rnemd) {
    print "UWAGA: liczba chunkow profilu (",nchunks,") != Nbin RNEMD (",Nbin_rnemd,")."
}
if (profiles_per_block < min_profiles_per_block) {
    print "UWAGA: blok ma tylko ok. ",profiles_per_block," profili T; zwieksz block_ps lub zapisuj profil czesciej."
}

step_of_row(r) = profile_first + int(r/(nchunks+1))*profile_every

if (auto_fit_range == 1) {
    if (nchunks >= 12) {
        fit_L_first_bin = 3
        fit_L_last_bin  = int(nchunks/2)-1
        fit_R_first_bin = int(nchunks/2)+3
        fit_R_last_bin  = nchunks-1
    } else {
        fit_L_first_bin = 2
        fit_L_last_bin  = int(nchunks/2)
        fit_R_first_bin = int(nchunks/2)+2
        fit_R_last_bin  = nchunks
    }
}

xLmin = (fit_L_first_bin-1)*dx
xLmax = fit_L_last_bin*dx
xRmin = (fit_R_first_bin-1)*dx
xRmax = fit_R_last_bin*dx

set table "ANALIZA_profil_sredni.dat"
plot profile using \
(($1>=1 && $1<=nchunks && step_of_row($0)>=avg_start && step_of_row($0)<=avg_end) ? $2*Lx : 1/0):\
(($1>=1 && $1<=nchunks && step_of_row($0)>=avg_start && step_of_row($0)<=avg_end) ? $4 : 1/0) \
smooth unique notitle
unset table

# regresja liniowa dwoch galezi
stats [xLmin:xLmax] "ANALIZA_profil_sredni.dat" using 1:2 name "LSTAT" nooutput
stats [xRmin:xRmax] "ANALIZA_profil_sredni.dat" using 1:2 name "RSTAT" nooutput

aL = LSTAT_slope
bL = LSTAT_intercept
aR = RSTAT_slope
bR = RSTAT_intercept

aL_error = LSTAT_slope_err
aR_error = RSTAT_slope_err

fL(x)=aL*x+bL
fR(x)=aR*x+bR

gradient_L = abs(aL)
gradient_R = abs(aR)
gradient   = (gradient_L+gradient_R)/2.0
gradient_fit_unc = 0.5*sqrt(aL_error**2+aR_error**2)
gradient_asym = 200.0*abs(gradient_L-gradient_R)/(gradient_L+gradient_R)

R2_L = LSTAT_correlation**2
R2_R = RSTAT_correlation**2
nfit_L = LSTAT_records
nfit_R = RSTAT_records

stats "ANALIZA_profil_sredni.dat" using 2 name "TPROF" nooutput
T_profile_min = TPROF_min
T_profile_max = TPROF_max
DeltaT_profile = T_profile_max-T_profile_min
DeltaT_rel_percent = 100.0*DeltaT_profile/rnemd_T_mean

# temperatury dwoch warstw wymiany
bin_1   = 1
bin_mid = int(nchunks/2)+1

stats profile using \
(($1==bin_1 && step_of_row($0)>=avg_start && step_of_row($0)<=avg_end) ? $4 : 1/0) \
name "TBIN1" nooutput

stats profile using \
(($1==bin_mid && step_of_row($0)>=avg_start && step_of_row($0)<=avg_end) ? $4 : 1/0) \
name "TBINM" nooutput

T_bin1_mean = TBIN1_mean
T_binmid_mean = TBINM_mean
DeltaT_exchange_mean = abs(T_binmid_mean-T_bin1_mean)

# STRUMIEN ENERGII I KAPPA

stats logfile using \
((valid(8) && $1>=avg_start && $1<=avg_end) ? ($1-avg_start)*h_ps : 1/0):\
((valid(8) && $1>=avg_start && $1<=avg_end) ? $5*heat_to_eV : 1/0) \
name "ESTAT" nooutput

mE       = ESTAT_slope
mE_error = ESTAT_slope_err
bE       = ESTAT_intercept
Efit(x)  = mE*x+bE

E_rate_fit = abs(mE)
R2_E = ESTAT_correlation**2

# kontrola DeltaE/t; poprawna dla monotonicznego skumulowanego f_heat
stats logfile using \
((valid(8) && $1>=avg_start && $1<=avg_end) ? $5*heat_to_eV : 1/0) \
name "HEAT" nooutput

DeltaE = HEAT_max-HEAT_min
E_rate_delta = DeltaE/t_avg_ps
E_rate_diff_percent = 100.0*abs(E_rate_fit-E_rate_delta)/E_rate_fit

# wzor Fouriera dla k
conv_k = 1602.176634
J_eV_ps_A2 = E_rate_fit/(2.0*A_A2)
J_W_m2 = J_eV_ps_A2*1.602176634e13

kappa_L = E_rate_fit/(2.0*A_A2*gradient_L)*conv_k
kappa_R = E_rate_fit/(2.0*A_A2*gradient_R)*conv_k
kappa_fit = E_rate_fit/(2.0*A_A2*gradient)*conv_k
kappa_fit_unc = kappa_fit*sqrt((mE_error/E_rate_fit)**2 + \
                               (gradient_fit_unc/gradient)**2)

# przewodnictwo arkuszowe G_2D = kappa*d
G2D_W_K  = kappa_fit*d*1.0e-10
G2D_nW_K = G2D_W_K*1.0e9

# BLOKI STEADY-STATE

nblocks = int(Nsteps_avg/block_steps)
unused_tail_ps = (Nsteps_avg-nblocks*block_steps)*h_ps

set print "ANALIZA_kappa_bloki.dat"
print "#bstart bend mid_step time_mid_ps grad_L grad_R grad_avg asym_percent kappa R2_L R2_R R2_E"

if (nblocks >= 1) {
    do for [ib=0:nblocks-1] {
        bstart = avg_start+ib*block_steps
        bend   = bstart+block_steps

        set table "ANALIZA_profil_tmp.dat"
        plot profile using \
        (($1>=1 && $1<=nchunks && step_of_row($0)>bstart && step_of_row($0)<=bend) ? $2*Lx : 1/0):\
        (($1>=1 && $1<=nchunks && step_of_row($0)>bstart && step_of_row($0)<=bend) ? $4 : 1/0) \
        smooth unique notitle
        unset table

        stats [xLmin:xLmax] "ANALIZA_profil_tmp.dat" using 1:2 name "BL" nooutput
        stats [xRmin:xRmax] "ANALIZA_profil_tmp.dat" using 1:2 name "BR" nooutput
        bgL = abs(BL_slope)
        bgR = abs(BR_slope)
        bg  = (bgL+bgR)/2.0
        basym = 200.0*abs(bgL-bgR)/(bgL+bgR)

        stats logfile using \
        ((valid(8) && $1>bstart && $1<=bend) ? ($1-rnemd_start)*h_ps : 1/0):\
        ((valid(8) && $1>bstart && $1<=bend) ? $5*heat_to_eV : 1/0) \
        name "BE" nooutput

        bErate = abs(BE_slope)
        bkappa = bErate/(2.0*A_A2*bg)*conv_k
        bmid = (bstart+bend)/2.0
        print bstart,bend,bmid,(bmid-rnemd_start)*h_ps,bgL,bgR,bg,basym,bkappa,BL_correlation**2,BR_correlation**2,BE_correlation**2
    }
}
set print

kappa_block_mean = kappa_fit
kappa_block_std  = 0.0
kappa_block_sem  = 0.0
kappa_SEM_percent = 999.0
kappa_SD_percent  = 999.0
kappa_block_agreement_percent = 999.0
kappa_block_trend_percent = 999.0

if (nblocks >= 1) {
    stats "ANALIZA_kappa_bloki.dat" using 9 name "KBLOCK" nooutput
    kappa_block_mean = KBLOCK_mean
    kappa_block_std  = (nblocks > 1) ? KBLOCK_ssd : 0.0
    kappa_block_sem  = (nblocks > 1) ? kappa_block_std/sqrt(nblocks) : 0.0
    kappa_block_agreement_percent = 100.0*abs(kappa_block_mean-kappa_fit)/abs(kappa_fit)

    if (nblocks > 1 && kappa_block_mean != 0) {
        kappa_SEM_percent = 100.0*kappa_block_sem/abs(kappa_block_mean)
        kappa_SD_percent  = 100.0*kappa_block_std/abs(kappa_block_mean)
        stats "ANALIZA_kappa_bloki.dat" using 4:9 name "KBTREND" nooutput
        block_span_ps = (nblocks-1)*block_ps
        kappa_block_trend_percent = 100.0*abs(KBTREND_slope)*block_span_ps/abs(kappa_block_mean)
    }
}

# KONWERGENCJA SKUMULOWANA

set print "ANALIZA_kappa_skumulowana.dat"
print "#kend time_ps grad_L grad_R grad_avg asym_percent kappa"

ncum = int(ceil(Nsteps_avg*1.0/conv_every))
do for [ic=1:ncum] {
    kend = avg_start + ic*conv_every
    if (kend > avg_end) {kend = avg_end}

    set table "ANALIZA_profil_tmp.dat"
    plot profile using \
    (($1>=1 && $1<=nchunks && step_of_row($0)>=avg_start && step_of_row($0)<=kend) ? $2*Lx : 1/0):\
    (($1>=1 && $1<=nchunks && step_of_row($0)>=avg_start && step_of_row($0)<=kend) ? $4 : 1/0) \
    smooth unique notitle
    unset table

    stats [xLmin:xLmax] "ANALIZA_profil_tmp.dat" using 1:2 name "CL" nooutput
    stats [xRmin:xRmax] "ANALIZA_profil_tmp.dat" using 1:2 name "CR" nooutput

    cgL = abs(CL_slope)
    cgR = abs(CR_slope)
    cg  = (cgL+cgR)/2.0
    casym = 200.0*abs(cgL-cgR)/(cgL+cgR)

    stats logfile using \
    ((valid(8) && $1>=avg_start && $1<=kend) ? ($1-rnemd_start)*h_ps : 1/0):\
    ((valid(8) && $1>=avg_start && $1<=kend) ? $5*heat_to_eV : 1/0) \
    name "CE" nooutput

    cErate = abs(CE_slope)
    ckappa = cErate/(2.0*A_A2*cg)*conv_k

    print kend,(kend-rnemd_start)*h_ps,cgL,cgR,cg,casym,ckappa
}
set print

cum_tail_start_ps = (avg_start-rnemd_start)*h_ps + 0.70*t_avg_ps
stats "ANALIZA_kappa_skumulowana.dat" using \
(($2>=cum_tail_start_ps) ? $2 : 1/0):\
(($2>=cum_tail_start_ps) ? $7 : 1/0) name "KCUM" nooutput

if (KCUM_records >= 2 && KCUM_mean_y != 0) {
    kappa_cum_trend_percent = 100.0*abs(KCUM_slope)*(0.30*t_avg_ps)/abs(KCUM_mean_y)
} else {
    kappa_cum_trend_percent = 999.0
}

# STATUS I WYNIK

heat_ok    = (R2_E >= heat_R2_min)
profile_ok = (R2_L >= profile_R2_min && R2_R >= profile_R2_min)
asym_ok    = (gradient_asym <= gradient_asym_max_percent)
blocks_ok  = (nblocks >= min_blocks_for_stats && profiles_per_block >= min_profiles_per_block)
sem_ok     = (blocks_ok && kappa_SEM_percent <= kappa_SEM_max_percent)
trend_ok   = (blocks_ok && kappa_block_trend_percent <= kappa_trend_max_percent && \
              kappa_cum_trend_percent <= kappa_trend_max_percent)
block_agreement_ok = (blocks_ok && kappa_block_agreement_percent <= kappa_block_agreement_max_percent)
deltaT_warn = (DeltaT_rel_percent > deltaT_rel_warn_percent)

rnemd_ok = (heat_ok && profile_ok && asym_ok && rnemd_temp_ok && energy_ok && \
            sem_ok && trend_ok && block_agreement_ok)

# glowna kappa = cale steady-state; SEM blokow jako niepewnosc statystyczna
kappa_final = kappa_fit
kappa_final_unc = (blocks_ok) ? kappa_block_sem : kappa_fit_unc
G2D_final_nW_K = kappa_final*d*0.1
G2D_unc_nW_K   = kappa_final_unc*d*0.1

if (kappa_SEM_percent <= 2.0) {
    sem_quality="BARDZO DOBRA"
} else {
    if (kappa_SEM_percent <= 5.0) {
        sem_quality="DOBRA"
    } else {
        if (kappa_SEM_percent <= 10.0) {sem_quality="UMIARKOWANA"}
        else {
            sem_quality="SLABA"
            }
    }
}

if (gradient_asym <= 5.0) {
    asym_quality="BARDZO DOBRA"
} else {
    if (gradient_asym <= 10.0) {asym_quality="AKCEPTOWALNA"}
    else {
        asym_quality="UWAGA"
        }
}

status_NVT="UWAGA"
status_NPT="UWAGA"
status_RNEMD="UWAGA"
if (nvt_ok) status_NVT="OK"
if (npt_ok) status_NPT="OK"
if (rnemd_ok) status_RNEMD="OK"

# temperatury warstw wymiany
set table $BIN1
plot profile using \
((step_of_row($0)-rnemd_start)*h_ps):4 \
with table if ($1==bin_1)

set table $BINMID
plot profile using \
((step_of_row($0)-rnemd_start)*h_ps):4 \
with table if ($1==bin_mid)

unset table

# WYKRESY

set grid back lc rgb "#D9DDE1"
set border linewidth 1.2 lc rgb "#333333"
set tics out nomirror
set samples 500
set key top center horizontal opaque box samplen 2 spacing 1.1

set style line 1 lc rgb "#6f2dbd" lw 2.0 pt 7 ps 0.9
set style line 2 lc rgb "#0f8a73" lw 2.2 pt 5 ps 0.9
set style line 3 lc rgb "#2F9ED8" lw 2.2 pt 7 ps 0.9
set style line 4 lc rgb "#D55E00" lw 2.0 dt 2
set style line 5 lc rgb "#555555" lw 1.6 dt 3

t_npt_boundary = (npt_start-nvt_start)*h_ps
t_avg_start = (avg_start-rnemd_start)*h_ps
t_avg_end   = (avg_end-rnemd_start)*h_ps

do for [fmt in "png svg"] {

    if (fmt eq "png") {
        set terminal pngcairo size 1600,900 enhanced font "Arial,13"
    } else {
        set terminal svg size 1600,900 enhanced font "Arial,13"
    }

    # equilibracja
    set output sprintf("%s_ANALIZA_01_equilibracja.%s",tag,fmt)
    set multiplot layout 2,1 title "Equilibracja NVT + NPT" font "Arial,16"

    set xlabel ""
    set ylabel "T [K]"
    set xrange [0:t_nvt_ps+t_npt_ps]
    set arrow 101 from t_npt_boundary,graph 0 to t_npt_boundary,graph 1 nohead dt 2 lc rgb "#777777"

    plot logfile using \
    ((valid(6) && $3>0 && $1>=nvt_start && $1<nvt_end) ? ($1-nvt_start)*h_ps : 1/0):5 \
    with lines ls 1 title "NVT", \
    logfile using \
    ((valid(6) && $3>0 && $1>npt_start && $1<npt_end) ? ($1-nvt_start)*h_ps : 1/0):5 \
    with lines ls 3 title "NPT", \
    T_target with lines ls 5 title sprintf("T zadana = %.1f K",T_target)

    set xlabel "czas od startu NVT [ps]"
    set ylabel "PE"

    plot logfile using \
    ((valid(6) && $3>0 && $1>=nvt_start && $1<nvt_end) ? ($1-nvt_start)*h_ps : 1/0):2 \
    with lines ls 1 title "PE NVT", \
    logfile using \
    ((valid(6) && $3>0 && $1>npt_start && $1<npt_end) ? ($1-nvt_start)*h_ps : 1/0):2 \
    with lines ls 3 title "PE NPT"

    unset arrow 101
    unset multiplot
    unset output
    set xrange [*:*]

    # NPT geometria
    if (npt_has_dimensions == 1) {
        set output sprintf("%s_ANALIZA_02_NPT_geometria.%s",tag,fmt)
        set title sprintf("NPT - geometria | dryf L_x = %.3f%%, L_y = %.3f%%",npt_Lx_drift_percent,npt_Ly_drift_percent)
        set xlabel "czas NPT [ps]"
        set ylabel "L [A]"

        plot logfile using \
        ((valid(10) && $3>0 && $1>npt_start && $1<npt_end) ? ($1-npt_start)*h_ps : 1/0):9 \
        with lines ls 2 title "L_x", \
        logfile using \
        ((valid(10) && $3>0 && $1>npt_start && $1<npt_end) ? ($1-npt_start)*h_ps : 1/0):10 \
        with lines ls 3 title "L_y"

        unset output
    }

    # RNEMD temperatura globalna
    set output sprintf("%s_ANALIZA_03_RNEMD_T.%s",tag,fmt)
    set title sprintf("RNEMD - T | srednia = %.2f K, dryf = %.3f%%",rnemd_T_mean,rnemd_T_drift_percent)
    set xlabel "czas RNEMD [ps]"
    set ylabel "T [K]"
    set xrange [0:t_rnemd_ps]
    set object 203 rect from first t_avg_start,graph 0 to first t_avg_end,graph 1 behind \
        fc rgb "#EEF4F7" fs solid 1.0 noborder

    plot logfile using \
    ((valid(8) && $1>rnemd_start && $1<=rnemd_end) ? ($1-rnemd_start)*h_ps : 1/0):3 \
    with lines ls 1 title "T RNEMD", \
    T_target with lines ls 5 title sprintf("T zadana = %.1f K",T_target), \
    rnemd_T_mean with lines ls 4 title sprintf("srednia steady-state = %.2f K",rnemd_T_mean)

    unset object 203
    unset output
    set xrange [*:*]

    # RNEMD temperatury warstw wymiany
    set output sprintf("%s_ANALIZA_04_RNEMD_warstwy.%s",tag,fmt)
    set title sprintf("RNEMD - warstwy wymiany | srednia Delta T = %.2f K",DeltaT_exchange_mean)
    set xlabel "czas RNEMD [ps]"
    set ylabel "T [K]"
    set xrange [0:t_rnemd_ps]


    plot $BIN1 using 1:2 with lines ls 1 title sprintf("T_{bin %d}",bin_1), \
     $BINMID using 1:2 with lines ls 3 title sprintf("T_{bin %d}",bin_mid)

    unset output
    set xrange [*:*]

    # f_heat
    set output sprintf("%s_ANALIZA_05_RNEMD_fheat.%s",tag,fmt)
    set title sprintf("RNEMD - f_{heat} | dE/dt = %.6f +/- %.6f eV/ps, R^2 = %.6f",E_rate_fit,mE_error,R2_E)
    set xlabel "czas RNEMD [ps]"
    set ylabel "f_{heat} [eV]"
    set xrange [0:t_rnemd_ps]

    plot logfile using \
    ((valid(8) && $1>rnemd_start && $1<=rnemd_end) ? ($1-rnemd_start)*h_ps : 1/0):($5*heat_to_eV) \
    with lines ls 2 title "f_{heat}", \
    [t_avg_start:t_avg_end] Efit(x-t_avg_start) with lines ls 4 title "fit steady-state"

    unset output
    set xrange [*:*]

    # gradient temperatury
    set output sprintf("%s_ANALIZA_06_gradient_T.%s",tag,fmt)
    set title sprintf("Profil T | |dT/dx| = %.5f K/A, asymetria = %.2f%%",gradient,gradient_asym)
    set xlabel "x [A]"
    set ylabel "T [K]"

    plot "ANALIZA_profil_sredni.dat" using 1:2 with points ls 1 title "T(x) - srednia", \
    [xLmin:xLmax] fL(x) with lines ls 2 \
    title sprintf("L: %.5f +/- %.5f K/A, R^2=%.4f",gradient_L,aL_error,R2_L), \
    [xRmin:xRmax] fR(x) with lines ls 3 \
    title sprintf("R: %.5f +/- %.5f K/A, R^2=%.4f",gradient_R,aR_error,R2_R)

    unset output

    # gradient w blokach
    set output sprintf("%s_ANALIZA_07_gradient_czas.%s",tag,fmt)
    set title sprintf("Gradient w blokach | blok = %.0f ps",block_ps)
    set xlabel "czas RNEMD [ps]"
    set ylabel "|dT/dx| [K/A]"
    set xrange [0:t_rnemd_ps]

    plot "ANALIZA_kappa_bloki.dat" using 4:5 with linespoints ls 1 title "gradient L", \
    "ANALIZA_kappa_bloki.dat" using 4:6 with linespoints ls 2 title "gradient R", \
    "ANALIZA_kappa_bloki.dat" using 4:7 with linespoints ls 3 title "gradient sredni"

    unset output
    set xrange [*:*]

    # kappa i konwergencja
    set output sprintf("%s_ANALIZA_08_kappa_czas.%s",tag,fmt)
    set title sprintf("Stabilnosc i konwergencja kappa | kappa = %.3f +/- %.3f W/(m K)",kappa_final,kappa_final_unc)
    set xlabel "czas RNEMD [ps]"
    set ylabel "kappa [W/(m K)]"
    set xrange [0:t_rnemd_ps]

    plot "ANALIZA_kappa_bloki.dat" using 4:9 with linespoints ls 1 title "kappa - bloki", \
    "ANALIZA_kappa_skumulowana.dat" using 2:7 with linespoints ls 2 title "kappa - skumulowana od avg_{start}", \
    [t_avg_start:t_avg_end] kappa_final with lines ls 4 title sprintf("wynik = %.3f +/- %.3f",kappa_final,kappa_final_unc)

    unset output
    set xrange [*:*]
    unset title
}

# WYNIKI TXT

wyniki_txt = sprintf("%s_ANALIZA_wyniki.txt",tag)
set print wyniki_txt

print "============================================================"
print "ANALIZA SYMULACJI RNEMD"
print "============================================================"
print ""

print "PARAMETRY"
print "material = ",material
print "potencjal = ",potential
print "anizotropia = ",anisotropy
print "T zadana [K] = ",T_target
print "timestep [ps] = ",h_ps
print "NVT [ps] = ",t_nvt_ps
print "NPT [ps] = ",t_npt_ps
print "RNEMD [ps] = ",t_rnemd_ps
print "steady-state kroki = ",avg_start," - ",avg_end
print "steady-state [ps] = ",t_avg_ps
print "swap co [krokow] = ",swap_every
print "swap co [ps] = ",swap_every*h_ps
print "Nbin RNEMD = ",Nbin_rnemd
print "plastry profilu = ",nchunks
print "dx [A] = ",dx
print "profile w steady-state = ",nprofiles_avg
print "blok statystyczny [ps] = ",block_ps
print "profile na blok ~= ",profiles_per_block
print ""

print "GEOMETRIA"
print "Lx [A] = ",Lx
print "Ly [A] = ",Ly
print "d [A] = ",d
print "A=Ly*d [A^2] = ",A_A2
print "SD Lx podczas RNEMD [A] = ",LXRN_stddev
print "SD Ly podczas RNEMD [A] = ",LYRN_stddev
if (npt_has_dimensions == 1) {
    print "Lx NPT srednia [A] = ",npt_Lx_mean
    print "Ly NPT srednia [A] = ",npt_Ly_mean
    print "dryf Lx NPT [%] = ",npt_Lx_drift_percent
    print "dryf Ly NPT [%] = ",npt_Ly_drift_percent
    print "Pxx NPT [atm] = ",npt_Pxx_mean," +/- ",npt_Pxx_std
    print "Pyy NPT [atm] = ",npt_Pyy_mean," +/- ",npt_Pyy_std
}
print ""

print "STABILNOSC NVT"
print "status = ",status_NVT
print "T [K] = ",nvt_T_mean," +/- ",nvt_T_std
print "odchylenie T od zadanej [%] = ",nvt_T_dev_percent
print "dryf T [% okna] = ",nvt_T_drift_percent
print "trend PE [liczba SD/okno] = ",nvt_PE_shift_sigma
print ""

print "STABILNOSC NPT"
print "status = ",status_NPT
print "T [K] = ",npt_T_mean," +/- ",npt_T_std
print "odchylenie T od zadanej [%] = ",npt_T_dev_percent
print "dryf T [% okna] = ",npt_T_drift_percent
print "Press [atm] = ",npt_P_mean," +/- ",npt_P_std
print ""

print "STABILNOSC RNEMD"
print "status = ",status_RNEMD
print "T [K] = ",rnemd_T_mean," +/- ",rnemd_T_std
print "odchylenie T od zadanej [%] = ",rnemd_T_dev_percent
print "dryf T [% okna] = ",rnemd_T_drift_percent
if (rnemd_has_etot == 1) print "trend Etot [liczba SD/okno] = ",rnemd_E_shift_sigma
print "T bin 1 [K] = ",T_bin1_mean
print "T bin srodkowy [K] = ",T_binmid_mean
print "DeltaT warstw wymiany [K] = ",DeltaT_exchange_mean
print ""

print "STRUMIEN ENERGII"
print "R2 f_heat = ",R2_E
print "dE/dt [eV/ps] = ",E_rate_fit
print "blad dE/dt [eV/ps] = ",mE_error
print "kontrola DeltaE/t [eV/ps] = ",E_rate_delta
print "roznica metod [%] = ",E_rate_diff_percent
print "J [W/m^2] = ",J_W_m2
print ""

print "GRADIENT TEMPERATURY"
print "fit L: plastry ",fit_L_first_bin," - ",fit_L_last_bin,"; punkty = ",nfit_L
print "fit R: plastry ",fit_R_first_bin," - ",fit_R_last_bin,"; punkty = ",nfit_R
print "|dT/dx| L [K/A] = ",gradient_L," +/- ",aL_error
print "R2 L = ",R2_L
print "|dT/dx| R [K/A] = ",gradient_R," +/- ",aR_error
print "R2 R = ",R2_R
print "sredni |dT/dx| [K/A] = ",gradient
print "asymetria L/R [%] = ",gradient_asym," (",asym_quality,")"
print "DeltaT profilu [K] = ",DeltaT_profile
print "DeltaT/T [%] = ",DeltaT_rel_percent
print ""

print "PRZEWODNICTWO CIEPLNE"
print "kappa L [W/(m K)] = ",kappa_L
print "kappa R [W/(m K)] = ",kappa_R
print "kappa steady-state [W/(m K)] = ",kappa_fit
print "niepewnosc tylko z fitow [W/(m K)] = ",kappa_fit_unc
print "niepokrywajace sie bloki w steady-state = ",nblocks
print "kappa srednia z blokow [W/(m K)] = ",kappa_block_mean
print "SD blokow [W/(m K)] = ",kappa_block_std
print "SEM blokow [W/(m K)] = ",kappa_block_sem
print "SEM [%] = ",kappa_SEM_percent," (",sem_quality,")"
print "roznica full/bloki [%] = ",kappa_block_agreement_percent
print "trend blokow [% zakresu] = ",kappa_block_trend_percent
print "trend konca kappa skumulowanej [%] = ",kappa_cum_trend_percent
print ""
print "WYNIK DO TABELI"
print "kappa [W/(m K)] = ",kappa_final," +/- ",kappa_final_unc
print "G_2D [nW/K] = ",G2D_final_nW_K," +/- ",G2D_unc_nW_K
print "UWAGA: blad nie obejmuje niepewnosci potencjalu, definicji d, finite-size ani metody RNEMD."
print ""

print "KRYTERIA QC"
print "NVT OK = ",nvt_ok
print "NPT termiczne OK = ",npt_thermal_ok
if (npt_has_dimensions == 1) print "NPT geometria OK = ",npt_geom_ok
print "RNEMD T OK = ",rnemd_temp_ok
if (rnemd_has_etot == 1) print "RNEMD Etot OK = ",energy_ok
print "f_heat liniowe = ",heat_ok
print "profil T liniowy L/R = ",profile_ok
print "asymetria gradientu OK = ",asym_ok
print "bloki statystyczne OK = ",blocks_ok
print "SEM OK = ",sem_ok
print "trend kappa OK = ",trend_ok
print "full/bloki zgodne = ",block_agreement_ok
print "duze DeltaT/T = ",deltaT_warn
print ""

print "INTERPRETACJA"
if (!nvt_ok) print "- NVT: koncowe okno nie spelnia ustawionych kryteriow."
if (!npt_thermal_ok) print "- NPT: temperatura nadal wykazuje odchylenie lub dryf."
if (npt_has_dimensions == 1 && !npt_geom_ok) print "- NPT: Lx/Ly nadal wykazuja trend."
if (!rnemd_temp_ok) print "- RNEMD: srednia T lub jej dryf przekracza kryterium."
if (rnemd_has_etot == 1 && !energy_ok) print "- RNEMD: Etot wykazuje trend; sprawdz timestep/integracje/QEq."
if (!heat_ok) print "- RNEMD: f_heat nie jest dostatecznie liniowe."
if (!profile_ok) print "- RNEMD: fit profilu T jest slaby; sprawdz zakres fitu i steady-state."
if (!asym_ok) print "- RNEMD: gradient L/R jest asymetryczny."
if (!blocks_ok) print "- RNEMD: za malo blokow lub profili na blok do mocnego SEM."
if (blocks_ok && !sem_ok) print "- RNEMD: SEM jest duzy."
if (blocks_ok && !trend_ok) print "- RNEMD: kappa nadal wykazuje trend."
if (blocks_ok && !block_agreement_ok) print "- RNEMD: kappa full i srednia blokowa roznia sie za mocno."
if (deltaT_warn) print "- RNEMD: DeltaT/T jest duze; warto przetestowac slabsze wymuszenie."
if (rnemd_ok) print "- RNEMD: wybrany zakres spelnia ustawione kryteria robocze."

set print

summary_csv = sprintf("%s_ANALIZA_summary.csv",tag)
set print summary_csv
print "MXene;Anizotropia;Potencjal;kappa[W/(mK)];unc_kappa[W/(mK)];X[nm];Y[nm];Z[A];T_target[K];T_mean[K];DeltaT[K];RNEMD_steps;NVT_steps;NPT_steps;dt[ps];Nbin;dx[A];swap_steps;Lx[nm];Ly[nm];G2D[nW/K];dEdt[eV/ps];gradient[K/A];J[W/m2];R2heat;R2L;R2R;asym[%]"
print sprintf("%s;%s;%s;%.10g;%.10g;%s;%s;%.10g;%.10g;%.10g;%.10g;%d;%d;%d;%.10g;%d;%.10g;%d;%.10g;%.10g;%.10g;%.10g;%.10g;%.10g;%.10g;%.10g;%.10g;%.10g", \
material,anisotropy,potential,kappa_final,kappa_final_unc,Lx_target_nm,Ly_target_nm,d,T_target,rnemd_T_mean,DeltaT_profile,Nsteps_rnemd,Nsteps_nvt,Nsteps_npt,h_ps,Nbin_rnemd,dx,swap_every,Lx_nm,Ly_nm,G2D_final_nW_K,E_rate_fit,gradient,J_W_m2,R2_E,R2_L,R2_R,gradient_asym)
set print

print ""
print "============================================================"
print "ANALIZA ZAKONCZONA"
print "Wyniki: ",wyniki_txt
print "CSV: ",summary_csv
print "============================================================"
