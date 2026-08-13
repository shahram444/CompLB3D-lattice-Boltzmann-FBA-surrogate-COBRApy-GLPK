/* ============================================================================
 * defineKinetics.hh  —  FULL ZHAO NETWORK : BIOTIC (MICROBIAL) REACTIONS
 *                        Option A  (12 biomass pools, reactions R74-R107)
 * ============================================================================
 *
 * HOW TO READ THIS FILE (for readers who do not program):
 *
 *   - This is a "recipe" file. It is not a program you run on its own. The main
 *     CompLaB3D solver opens it and calls the one routine inside (defineRxnKinetics)
 *     ONE time for every little cube of the domain (each cube is a "voxel"), on
 *     every time step. The solver hands the routine the local amounts of microbes
 *     and chemicals in that voxel, and the routine hands back how fast each of
 *     them is changing there.
 *
 *   - Anything after "//", or between the block-comment markers, is a note for
 *     humans. The computer ignores it. Read it as the explanation of the line.
 *
 *   - Numbers are passed as "lists" read by position, counting from 0. This file
 *     uses the 14 reacting chemicals C0..C13 and the 12 microbe pools B0..B11:
 *
 *        CHEMICALS  C0=acetate  C1=U(VI)   C2=U(IV)    C3=Fe(II)
 *                   C4=Fe(III) easy   C5=Fe(III) medium   C6=Fe(III) hard
 *                   C7=sulfate  C8=sulfide  C9=FeS(mackinawite)
 *                   C10=bicarbonate  C11=H+ (proton)  C12=methane  C13=S0 (elemental sulfur)
 *        (The solver's chemical vector also carries 10 equilibrium-only components
 *         C14..C23 - major ions and mineral surface sites - which this file never
 *         touches. This biotic file also leaves C13=S0 at 0; the abiotic file makes S0.)
 *
 *        MICROBES   B0=Geobacter/hard   B1=Geobacter/medium  B2=Geobacter/easy
 *                   B3=Geobacter/suspended
 *                   B4=Rhodoferax/hard  B5=Rhodoferax/medium B6=Rhodoferax/easy
 *                   B7=Rhodoferax/suspended
 *                   B8=methanogen/attached   B9=methanogen/suspended
 *                   B10=sulfate-reducer/attached  B11=sulfate-reducer/suspended
 *
 *   - Fe(III) oxide is split into three "crystallinity" classes because bacteria
 *     reduce them at different rates: easy (amorphous + clay-bound), medium
 *     (poorly crystalline), hard (goethite). This split is straight from Zhao's
 *     Supplementary Information (Table S1, footnotes c-e).
 *
 *   WHAT REACTIONS LIVE HERE (word form; solids marked (s)):
 *     Iron reducers grow by "breathing" solid Fe(III), releasing Fe(II):
 *        R74-76 (Geobacter)  acetate + Fe(III) -> biomass + Fe(II) + bicarbonate
 *        R93-95 (Rhodoferax) acetate + Fe(III) -> biomass + Fe(II) + bicarbonate
 *     Geobacter also turns dissolved U(VI) into trapped solid U(IV):
 *        R77-80  acetate + U(VI) -> biomass + U(IV)(s) + bicarbonate
 *        (this rate is throttled as FeS coats the surfaces - Zhao Eq. 1)
 *     Sulfate reducers grow on sulfate, releasing sulfide:
 *        R103-104 acetate + sulfate -> biomass + sulfide + bicarbonate
 *     Methanogens grow on acetate, making methane (inhibited by easy-Fe & sulfate):
 *        R88-89  acetate -> biomass + methane + CO2
 *     Plus first-order decay (R81-84, R91-92, R96-99, R105-106) and
 *     attachment/detachment between attached and suspended pools
 *     (R85-87, R90, R100-102, R107).
 *
 *   The abiotic (non-microbial) reactions R108-R114 live in the companion file
 *   defineAbioticKinetics.hh.  FeS (C9) is produced there, not here.
 *
 *   SAFEGUARDS CARRIED OVER FROM TEST T9 (see comments at each site):
 *     * every amount is floored at zero before use (no negative concentrations),
 *     * a single per-step "stability limit" scales ALL reactions by one common
 *       factor so no chemical is drawn below zero in a step AND element balances
 *       stay exact,
 *     * a final not-a-number / infinity guard zeroes any bad value,
 *     * dt_kinetics must equal the run's real ADE time step (see the note there).
 * ============================================================================ */
#ifndef DEFINE_KINETICS_HH   // include guard: the three #ifndef/#define/#endif lines
#define DEFINE_KINETICS_HH
/* ===========================================================================
 * [TRUE-UNITS 2026-07-24] ALL RATE CONSTANTS CONVERTED TO PHYSICAL PER-SECOND
 * (original Zhao/literature per-DAY value / 86400). Vmax: mol Ac.(mol cells)^-1.s^-1 ;
 * kd & attachment/detachment: s^-1. Ks (mol/L), Ki (mol/L), yields (mol/mol) and
 * stoichiometric coefficients are NOT time units and are unchanged.
 * dt_kinetics set to the run ade_dt (~3.0e-4 s). See parameter reference doc.
 * ============================================================================ */

#include <vector>       // the "list of numbers" type used to pass amounts in and rates out
#include <cmath>        // math helpers (for example, the not-a-number test used at the end)
#include <algorithm>    // small helpers such as std::max ("take the larger of two numbers")
#include <iostream>     // lets us print a warning to the screen if something goes wrong

/* ---------------------------------------------------------------------------
 * CHEMICAL and MICROBE POSITIONS (names for the list positions above)
 *   Using names instead of bare numbers makes the reactions readable and stops
 *   index mistakes. "constexpr int NAME = k;" just gives the number k a NAME.
 * --------------------------------------------------------------------------- */
namespace ZhaoIdx {
    // chemicals (positions in list C)
    constexpr int AC=0, U6=1, U4=2, FE2=3, FE3E=4, FE3M=5, FE3H=6,
                  SO4=7, HS=8, FES=9, HCO3=10, HP=11, CH4=12, S0=13;
    // microbes (positions in list B)
    constexpr int GEO_H=0, GEO_M=1, GEO_E=2, GEO_S=3,
                  RF_H=4,  RF_M=5,  RF_E=6,  RF_S=7,
                  MT_A=8,  MT_S=9,  SRB_A=10, SRB_S=11;
}

/* ---------------------------------------------------------------------------
 * THE FIXED NUMBERS (reaction constants)
 *   Rates and half-speed points from Zhao et al. Supplementary Information,
 *   Tables S3 and S4. Rate constants are PER-SECOND (Zhao per-day / 86400). The half-
 *   saturation constants Ks are in mol/L. "constexpr double NAME = VALUE;"
 *   means: a fixed decimal number we give a NAME so we can use it below.
 * --------------------------------------------------------------------------- */
namespace KineticParams {
    // -- Geobacter maximum specific acetate-uptake rates on each Fe class (Table S3 R74-R76)
    constexpr double Vmax_Geo_Feh = 8.68e-6;   // on hard (goethite) Fe(III)
    constexpr double Vmax_Geo_Fem = 1.74e-5;   // on medium (poorly crystalline) Fe(III)
    constexpr double Vmax_Geo_Fee = 1.74e-5;   // on easy (amorphous + clay) Fe(III)
    constexpr double Vmax_Geo_U   = 4.05e-6;   // Geobacter on U(VI)              (Table S3 R77-R80)
    // -- Rhodoferax maximum specific rate on Fe(III) (same for all classes)  (Table S3 R93-R95)
    constexpr double Vmax_Rf_Fe   = 2.08e-6;
    // -- methanogen and sulfate reducer maximum specific rates on acetate/sulfate
    constexpr double Vmax_Mt      = 1.62e-4;   // methanogen on acetate           (Table S3 R88-R89)
    constexpr double Vmax_SRB     = 2.55e-5;   // sulfate reducer on sulfate      (Table S3 R103-R104)

    // -- half-saturation constants (mol/L): the amount at which a rate is half its top speed
    constexpr double Ks_Fe        = 2.5e-3; // Fe(III) half-sat (Table S3 R75, medium; used for all classes)
    constexpr double Ks_U         = 2.0e-5; // U(VI) half-sat  (Table S3 R77-R80; Roden 2002, Istok 2004)
    constexpr double Ks_Ac_Geo    = 1.0e-4; // acetate half-sat, Geobacter (Table S3 R74-R80; Scheibe 2006)
    constexpr double Ks_Fe_Rf     = 2.5e-4; // Fe(III) half-sat, Rhodoferax (Table S3 R94; Barlett 2012)
    constexpr double Ks_Ac_Rf     = 5.4e-6; // acetate half-sat, Rhodoferax (Table S3 R93-R95; Zhuang 2011)
    constexpr double Ks_Ac_Mt     = 1.0e-5; // acetate half-sat, methanogen (Table S3 R88; Fukuzaki 1990)
    constexpr double Ks_SO4       = 2.0e-4; // sulfate half-sat, SRB (Table S3 R103; Ingvorsen 1984)
    constexpr double Ks_Ac_SRB    = 7.0e-5; // acetate half-sat, SRB (Ingvorsen 1984, acetate Km=70uM; CORRECTED from 2.0e-4 which was the sulfate Km)

    // -- methanogen inhibition constants (mol/L): methanogenesis is suppressed while
    //    easy Fe(III) or sulfate are present (Table S3 R88; Liu 2011, Oremland 1982)
    constexpr double Ki_Fee       = 1.0e-3; // inhibition by easy Fe(III)
    constexpr double Ki_SO4       = 1.6e-3; // inhibition by sulfate

    // -- yields Y (mol cells per mol acetate) (Table S4)
    constexpr double Y_Geo_Fe = 0.035;  // Geobacter on Fe(III) (Liu 2011, Scheibe 2009; midpoint 0.031-0.040)
    constexpr double Y_Geo_U  = 0.012;  // Geobacter on U(VI)   (Sanford 2007)
    constexpr double Y_Rf     = 0.062;  // Rhodoferax           (Zhuang 2011)
    constexpr double Y_Mt     = 0.034;  // methanogen           (Marchaim 1992)
    constexpr double Y_SRB    = 0.041;  // sulfate reducer      (Widdel & Pfenning 1981)

    // -- decay (die-off) rates kd (per second) (Table S3 R81-R84, R91-R92, R96-R99, R105-R106)
    constexpr double kd        = 6.94e-8; // most pools
    constexpr double kd_Geo_e  = 2.31e-8; // easy-Fe Geobacter decays slower (cells aggregate; England 1993)

    // -- attachment / detachment rates (per second) (Table S3 R85-R87, R90, R100-R102, R107)
    constexpr double k_att     = 1.62e-7;   // attachment (all; Zhao 2011)
    constexpr double k_det_h   = 2.08e-7;   // Geobacter detachment, hard  (Levy 2007)
    constexpr double k_det_m   = 2.55e-7;   // Geobacter detachment, medium
    constexpr double k_det_e   = 1.16e-8;   // Geobacter detachment, easy (aggregates; hard to leave)
    constexpr double k_det_Rf  = 1.62e-7;   // Rhodoferax detachment (h,m); easy uses k_det_e
    constexpr double k_det_Mt  = 6.94e-8;   // methanogen detachment
    constexpr double k_det_SRB = 1.62e-7;   // sulfate reducer detachment

    // -- Zhao stoichiometric coefficients per mol acetate (Supplementary Table S1) --
    constexpr double s_Fe_per_Ac   = 4.8;    // mol Fe(III) reduced per mol acetate (R74-76, R93-95)
    constexpr double s_HCO3_Geo_Fe = 1.845;  // bicarbonate released, Geobacter-Fe (R74-76)
    constexpr double s_Hp_Geo_Fe   = 6.17;   // protons consumed, Geobacter-Fe   (R74-76)
    constexpr double s_HCO3_Rf_Fe  = 1.69;   // bicarbonate released, Rhodoferax (R93-95)
    constexpr double s_Hp_Rf_Fe    = 6.94;   // protons consumed, Rhodoferax     (R93-95)
    constexpr double s_U_per_Ac    = 3.1;    // mol U(VI) reduced per mol acetate (R77-80)
    constexpr double s_HCO3_Geo_U  = 1.94;   // bicarbonate released, Geobacter-U (R77-80)
    constexpr double s_Hp_Geo_U    = 8.67;   // protons produced, Geobacter-U    (R77-80)
    constexpr double s_SO4_per_Ac  = 0.863;  // mol sulfate reduced per mol acetate (R103-104)
    constexpr double s_HCO3_SRB    = 1.80;   // bicarbonate released, SRB (R103-104)
    constexpr double s_Hp_SRB      = 0.22;   // protons consumed, SRB    (R103-104)
    constexpr double s_CH4_per_Ac  = 0.93;   // mol methane per mol acetate (R88-89)
    constexpr double s_CO2_per_Ac  = 0.93;   // mol CO2 (-> bicarbonate) per mol acetate (R88-89)
    constexpr double s_Hp_Mt       = 0.77;   // protons consumed, methanogen (R88-89)

    // -- housekeeping numbers --
    constexpr double eps        = 1.0e-12;  // tiny number so we never divide by zero
    constexpr double MIN_BIOMASS = 0.1;     // a voxel counts as "active" (reporting) above this biomass
    constexpr double MAX_RATE_FRACTION = 0.5;  // stability limit: a reaction may use at most half of a
                                               // chemical in one step (T9 safeguard)
    // IMPORTANT - TIME UNITS (matches the actual CompLaB3D solver). The solver advances the rates this
    // file returns using its reaction sub-step "ade_dt", which it computes IN SECONDS as
    //   ade_dt = ((tau-0.5)/3) * dx^2 / D0 ,   D0 = substrate-0 pore diffusion (acetate, 1e-9 m^2/s here).
    // For tau=0.8:  dx=1 um -> ade_dt = 1.0e-4 s ;  dx=2 um -> ade_dt = 4.0e-4 s. The Zhao rate constants
    // have been CONVERTED to true per-second (Zhao per-day / 86400); dt_kinetics equals the run ade_dt so the stability cap is physical.
    // dt_kinetics below is used ONLY by the stability limit and MUST equal the run's ade_dt (see the
    // "[ADE] dt" value printed in the log). Setting it too LARGE over-throttles growth (see the T9 notes);
    // using the true ade_dt makes the limit bind only when a step would genuinely over-draw a cell.
    constexpr double dt_kinetics = 1.0e-4;     // [seconds] = the run's printed [ADE] dt for this dx=1um duct (~1.0e-4). Set equal to ade_dt for any geometry.
}

/* ---------------------------------------------------------------------------
 * PROGRESS COUNTERS  (reporting only; no effect on the chemistry)
 * --------------------------------------------------------------------------- */
namespace KineticsStats {
    static double iter_sum_dB = 0.0, iter_max_biomass = 0.0, iter_max_dB = 0.0, iter_min_DOC = 1e30;
    static long   iter_cells_with_biomass = 0, iter_cells_with_growth = 0;
    inline void resetIteration() {                 // wipe the tallies clean at the start of each step
        iter_sum_dB=0; iter_max_biomass=0; iter_max_dB=0; iter_min_DOC=1e30;
        iter_cells_with_biomass=0; iter_cells_with_growth=0;
    }
    inline void accumulate(double biomass, double donor, double dB) {   // add one voxel to the tallies
        if (biomass > KineticParams::MIN_BIOMASS) {
            iter_cells_with_biomass++; iter_sum_dB += dB;
            if (biomass > iter_max_biomass) iter_max_biomass = biomass;
            if (dB > iter_max_dB) iter_max_dB = dB;
            if (donor < iter_min_DOC && donor > 0) iter_min_DOC = donor;
            if (dB > 0) iter_cells_with_growth++;
        }
    }
    inline void getStats(long& cb, long& cg, double& s, double& mB, double& mdB, double& mD) {
        cb=iter_cells_with_biomass; cg=iter_cells_with_growth;
        s=iter_sum_dB; mB=iter_max_biomass; mdB=iter_max_dB;
        mD=(iter_min_DOC<1e20)?iter_min_DOC:0.0;
    }
}

/* ---------------------------------------------------------------------------
 * SMALL HELPER: the Monod "throttle"  S/(Ks+S)
 *   A number that grows from 0 (no substrate) toward 1 (plenty of substrate).
 *   "inline double name(args){...}" defines a tiny reusable calculation.
 * --------------------------------------------------------------------------- */
inline double monod(double S, double Ks) { return S / (Ks + S); }

/* ---------------------------------------------------------------------------
 * THE MAIN REACTION ROUTINE
 *   Runs once per voxel per step. Inputs: B (12 microbe amounts) and C (the
 *   chemical amounts) here. Outputs (the "&" means the solver reads them back):
 *   subsR (how fast each chemical changes) and bioR (how fast each of the 12
 *   microbes changes). This file reads/writes only the 14 reacting chemicals C0..C13.
 * --------------------------------------------------------------------------- */
void defineRxnKinetics(
    std::vector<double> B,          // IN : 12 microbe amounts here (see the guide above)
    std::vector<double> C,          // IN : the chemical amounts here (reacting species C0..C13; see the guide above)
    std::vector<double>& subsR,     // OUT: fill with how fast each chemical changes (this file writes C0..C13)
    std::vector<double>& bioR,      // OUT: fill with how fast each of the 12 microbes changes
    plb::plint mask                 // IN : voxel label (pore/biofilm/solid); reactions only fire in pore/biofilm
) {
    using namespace KineticParams;
    using namespace ZhaoIdx;

    for (size_t i = 0; i < subsR.size(); ++i) subsR[i] = 0.0;   // "for" = repeat: start every chemical rate at 0
    for (size_t i = 0; i < bioR.size();  ++i) bioR[i]  = 0.0;   // start every microbe rate at 0
    if (B.size() < 12 || C.size() < 14) return;                 // if a list is short, "return" = stop, do nothing
    // note: this biotic file never makes elemental sulfur S0 (C13); the abiotic file does. subsR[S0] stays 0 here.
    if (mask < 2) return;                                       // do nothing in solid/bounce-back voxels

    // --- read the chemicals; std::max(x,0.0) floors any tiny negative at zero (T9 safeguard) ---
    double Ac   = std::max(C[AC],   0.0);   // acetate (the shared food / electron donor)
    double cU6  = std::max(C[U6],   0.0);   // dissolved U(VI)  (local name differs from the index U6)
    double Fe3e = std::max(C[FE3E], 0.0);   // solid Fe(III), easy class
    double Fe3m = std::max(C[FE3M], 0.0);   // solid Fe(III), medium class
    double Fe3h = std::max(C[FE3H], 0.0);   // solid Fe(III), hard class
    double SO4c = std::max(C[SO4],  0.0);   // sulfate (free SO4-- bucket for now; totalled below)
    double FeS  = std::max(C[FES],  0.0);   // FeS mineral formed so far (made in the abiotic file)

    // === FULL-NETWORK TOTALS (kinetics_on_totals; NO change needed anywhere else) ======
    //  With the full equilibrium speciation ON, C[AC]/C[U6]/C[SO4] hold only the FREE ion
    //  and the rest of the element sits in the complex buckets C[24..94]. Zhao's rate laws
    //  use the TOTAL DISSOLVED concentration, so rebuild it here by adding the AQUEOUS
    //  complexes (the immobile SURFACE-bound complexes are sorbed, not dissolved, so they
    //  are intentionally excluded). Indices are fixed by the 95-substrate full-network
    //  layout (species24..94, same order as CompLaB.xml; see Zhao Table S1). Guarded by
    //  C.size()>24 so this SAME file also runs unchanged with equilibrium OFF (24 subs).
    if (C.size() > 24) {
        Ac   += C[56]+C[81]+C[82]+C[83]+C[84]+C[85]+C[86];              // + NH4Ac, CH3COO-, Ca/Fe/Mg/Na/K-acetate
        cU6  += 2*(C[65]+C[66]+C[67]) + 3*(C[68]+C[69]+C[70]+C[71])      // + dimeric(x2)/trimeric(x3) U hydroxo-carbonates
              + C[72]+C[73]+C[74]+C[75]+C[76]+C[77]+C[78]+C[79]+C[80];   // + monomeric UO2 carbonate/hydroxo/phosphate
        SO4c += C[47]+C[48]+C[57];                                      // + KHSO4, KSO4-, NH4SO4-
    }
    // ================================================================================
    double Fe3tot = Fe3e + Fe3m + Fe3h;     // total Fe(III) surface still available

    // --- read the biomass pools; also floored at zero ---
    double Bgh=std::max(B[GEO_H],0.0), Bgm=std::max(B[GEO_M],0.0), Bge=std::max(B[GEO_E],0.0), Bgs=std::max(B[GEO_S],0.0);
    double Brh=std::max(B[RF_H],0.0),  Brm=std::max(B[RF_M],0.0),  Bre=std::max(B[RF_E],0.0),  Brs=std::max(B[RF_S],0.0);
    double Bma=std::max(B[MT_A],0.0),  Bms=std::max(B[MT_S],0.0);
    double Bsa=std::max(B[SRB_A],0.0), Bss=std::max(B[SRB_S],0.0);

    // --- shared "throttles" ---
    double m_ac_G  = monod(Ac, Ks_Ac_Geo);  // acetate limitation for Geobacter
    double m_ac_R  = monod(Ac, Ks_Ac_Rf);   // acetate limitation for Rhodoferax
    double m_ac_M  = monod(Ac, Ks_Ac_Mt);   // acetate limitation for methanogen
    double m_ac_S  = monod(Ac, Ks_Ac_SRB);  // acetate limitation for sulfate reducer
    double m_U     = monod(cU6, Ks_U);       // U(VI) limitation

    // Zhao Eq. 1 surface-competition factor: the share of cell surface still FREE of FeS.
    // ~1 when little FeS is present (mostly biological U reduction), ~0 once FeS coats the
    // surfaces (biological U reduction switches off; the abiotic file takes over via frac_ab).
    double frac_bio = 1.0 - FeS / (FeS + Fe3tot + eps);
    frac_bio = std::min(std::max(frac_bio, 0.0), 1.0);           // keep safely inside 0..1

    // methanogen inhibition: near 1 when easy-Fe and sulfate are low, near 0 when they are high
    double inhib_Mt = (Ki_Fee/(Ki_Fee + Fe3e)) * (Ki_SO4/(Ki_SO4 + SO4c));

    /* =====================================================================
     * REACTION RATES  (rho = specific acetate-uptake rate for that pathway,
     *   mol acetate per L per second = Vmax * biomass * throttles).  Chemical
     *   changes are then rho times Zhao's per-acetate stoichiometry, so every
     *   element closes by construction.
     * ===================================================================== */
    // -- Geobacter iron reduction on each Fe class (R74-76): acetate + Fe(III) -> Fe(II) + bicarbonate --
    double rGh_Fe = Vmax_Geo_Feh * Bgh * monod(Fe3h, Ks_Fe) * m_ac_G;   // on hard Fe(III)
    double rGm_Fe = Vmax_Geo_Fem * Bgm * monod(Fe3m, Ks_Fe) * m_ac_G;   // on medium Fe(III)
    double rGe_Fe = Vmax_Geo_Fee * Bge * monod(Fe3e, Ks_Fe) * m_ac_G;   // on easy Fe(III)

    // -- Geobacter U(VI) reduction (R77-80), throttled by the free-surface factor frac_bio (Zhao Eq.1) --
    double rGh_U = frac_bio * Vmax_Geo_U * Bgh * m_U * m_ac_G;          // hard-attached Geobacter on U(VI)
    double rGm_U = frac_bio * Vmax_Geo_U * Bgm * m_U * m_ac_G;          // medium-attached
    double rGe_U = frac_bio * Vmax_Geo_U * Bge * m_U * m_ac_G;          // easy-attached
    double rGs_U = frac_bio * Vmax_Geo_U * Bgs * m_U * m_ac_G;          // suspended Geobacter (R80)

    // -- Rhodoferax iron reduction on each Fe class (R93-95) --
    double rRh_Fe = Vmax_Rf_Fe * Brh * monod(Fe3h, Ks_Fe_Rf) * m_ac_R; // on hard
    double rRm_Fe = Vmax_Rf_Fe * Brm * monod(Fe3m, Ks_Fe_Rf) * m_ac_R; // on medium
    double rRe_Fe = Vmax_Rf_Fe * Bre * monod(Fe3e, Ks_Fe_Rf) * m_ac_R; // on easy

    // -- methanogen growth on acetate (R88-89), inhibited by easy-Fe and sulfate --
    double rMa = Vmax_Mt * Bma * m_ac_M * inhib_Mt;                     // attached methanogen
    double rMs = Vmax_Mt * Bms * m_ac_M * inhib_Mt;                     // suspended methanogen

    // -- sulfate reducer growth on sulfate (R103-104): acetate + sulfate -> sulfide + bicarbonate --
    double rSa = Vmax_SRB * Bsa * monod(SO4c, Ks_SO4) * m_ac_S;         // attached SRB
    double rSs = Vmax_SRB * Bss * monod(SO4c, Ks_SO4) * m_ac_S;         // suspended SRB

    /* ---------------------------------------------------------------------
     * PER-STEP STABILITY LIMIT (T9 safeguard)
     *   Acetate is shared by every pathway. Work out one shrink factor "s" so
     *   no chemical (acetate, each Fe class, U6, sulfate) loses more than half
     *   of what it holds in a step, then scale ALL reactions by that same "s".
     *   Scaling everything together keeps the element balances exact.
     * --------------------------------------------------------------------- */
    double s = 1.0;                                                     // 1.0 = no shrink needed yet
    double acDraw = rGh_Fe+rGm_Fe+rGe_Fe + rGh_U+rGm_U+rGe_U+rGs_U      // total acetate demanded this step
                  + rRh_Fe+rRm_Fe+rRe_Fe + rMa+rMs + rSa+rSs;
    double capAc  = Ac   * MAX_RATE_FRACTION / dt_kinetics;             // most acetate we allow to be used
    if (acDraw > capAc && acDraw > 0.0) s = std::min(s, capAc/acDraw);
    double fehDraw = s_Fe_per_Ac*(rGh_Fe + rRh_Fe);                     // hard Fe(III) consumed
    double capFeh  = Fe3h * MAX_RATE_FRACTION / dt_kinetics;
    if (fehDraw > capFeh && fehDraw > 0.0) s = std::min(s, capFeh/fehDraw);
    double femDraw = s_Fe_per_Ac*(rGm_Fe + rRm_Fe);                     // medium Fe(III) consumed
    double capFem  = Fe3m * MAX_RATE_FRACTION / dt_kinetics;
    if (femDraw > capFem && femDraw > 0.0) s = std::min(s, capFem/femDraw);
    double feeDraw = s_Fe_per_Ac*(rGe_Fe + rRe_Fe);                     // easy Fe(III) consumed
    double capFee  = Fe3e * MAX_RATE_FRACTION / dt_kinetics;
    if (feeDraw > capFee && feeDraw > 0.0) s = std::min(s, capFee/feeDraw);
    double uDraw   = s_U_per_Ac*(rGh_U+rGm_U+rGe_U+rGs_U);              // U(VI) consumed
    double capU    = cU6  * MAX_RATE_FRACTION / dt_kinetics;
    if (uDraw > capU && uDraw > 0.0) s = std::min(s, capU/uDraw);
    double soDraw  = s_SO4_per_Ac*(rSa+rSs);                            // sulfate consumed
    double capSO4  = SO4c * MAX_RATE_FRACTION / dt_kinetics;
    if (soDraw > capSO4 && soDraw > 0.0) s = std::min(s, capSO4/soDraw);

    // apply the one shrink factor to every reaction ("*= s" means multiply by s)
    rGh_Fe*=s; rGm_Fe*=s; rGe_Fe*=s;
    rGh_U*=s;  rGm_U*=s;  rGe_U*=s;  rGs_U*=s;
    rRh_Fe*=s; rRm_Fe*=s; rRe_Fe*=s;
    rMa*=s; rMs*=s; rSa*=s; rSs*=s;

    // convenient sums
    double GeoFe = rGh_Fe + rGm_Fe + rGe_Fe;          // all Geobacter iron reduction
    double RfFe  = rRh_Fe + rRm_Fe + rRe_Fe;          // all Rhodoferax iron reduction
    double GeoU  = rGh_U + rGm_U + rGe_U + rGs_U;      // all biotic U(VI) reduction
    double Mt    = rMa + rMs;                          // all methanogenesis
    double SRB   = rSa + rSs;                          // all sulfate reduction

    /* ---------------------------------------------------------------------
     * WRITE HOW FAST EACH CHEMICAL CHANGES (rate = rho * per-acetate stoichiometry)
     * --------------------------------------------------------------------- */
    // helper: remove a total consumption 'cons' (mol/L/s) of a dissolved element by
    // spreading it over its buckets IN PROPORTION TO CONTENT, so the element is conserved
    // and no bucket is driven negative even when the free ion is a tiny fraction of the
    // total. subsR[i] += -cons*C[i]/tot for every listed bucket i; sum over buckets
    // (weighted by their stoichiometry) removes exactly 'cons'. With equilibrium OFF the
    // list collapses to the primary bucket and this equals the old single-bucket write.
    auto drawFrom = [&](double cons, double tot, std::initializer_list<int> idx){
        if (tot <= 0.0 || cons <= 0.0) return;
        for (int i : idx) if (i < (int)C.size()) subsR[i] += -cons * std::max(C[i],0.0) / tot;
    };

    // acetate: consumed by every growth pathway (one mol acetate per unit rho)
    double acCons = GeoFe + GeoU + RfFe + Mt + SRB;
    drawFrom(acCons, Ac, {AC,56,81,82,83,84,85,86});
    // iron: Fe(III) of each class consumed by iron reducers; Fe(II) released (4.8 per acetate)
    subsR[FE3H] = -s_Fe_per_Ac * (rGh_Fe + rRh_Fe);
    subsR[FE3M] = -s_Fe_per_Ac * (rGm_Fe + rRm_Fe);
    subsR[FE3E] = -s_Fe_per_Ac * (rGe_Fe + rRe_Fe);
    subsR[FE2]  = +s_Fe_per_Ac * (GeoFe + RfFe);       // dissolved Fe(II) produced (abiotic file consumes it)
    // uranium: U(VI) consumed (spread over free + aqueous U complexes), immobile U(IV) produced (3.1 per acetate)
    double uCons = s_U_per_Ac * GeoU;
    drawFrom(uCons, cU6, {U6,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80});
    subsR[U4]  = +s_U_per_Ac * GeoU;
    // sulfur: sulfate consumed (spread over free + sulfate complexes), sulfide released (0.863 per acetate)
    double soCons = s_SO4_per_Ac * SRB;
    drawFrom(soCons, SO4c, {SO4,47,48,57});
    subsR[HS]  = +s_SO4_per_Ac * SRB;
    // methane produced by methanogens (0.93 per acetate)
    subsR[CH4] = +s_CH4_per_Ac * Mt;
    // bicarbonate released by every pathway (its own per-acetate coefficient)
    subsR[HCO3] = s_HCO3_Geo_Fe*GeoFe + s_HCO3_Geo_U*GeoU + s_HCO3_Rf_Fe*RfFe
                + s_CO2_per_Ac*Mt   + s_HCO3_SRB*SRB;
    // protons: consumed by iron/sulfate/methane pathways, produced by U reduction
    subsR[HP]  = -s_Hp_Geo_Fe*GeoFe - s_Hp_Rf_Fe*RfFe - s_Hp_SRB*SRB - s_Hp_Mt*Mt
                + s_Hp_Geo_U*GeoU;
    // subsR[FES] left at 0 here on purpose: FeS is made/consumed only in the abiotic file.

    /* ---------------------------------------------------------------------
     * WRITE HOW FAST EACH MICROBE CHANGES
     *   growth = yield * (acetate uptake of that pool), minus first-order decay,
     *   plus attachment/detachment transfers (which move biomass between the
     *   attached and suspended pools without creating or destroying any).
     * --------------------------------------------------------------------- */
    // attachment/detachment transfers (biomass-conserving): attachment from the single
    // suspended pool is split evenly across the three attached Fe-class pools.
    double gAtt_h = k_att*Bgs/3.0 - k_det_h*Bgh;   // Geobacter hard  <-> suspended
    double gAtt_m = k_att*Bgs/3.0 - k_det_m*Bgm;   // Geobacter medium<-> suspended
    double gAtt_e = k_att*Bgs/3.0 - k_det_e*Bge;   // Geobacter easy  <-> suspended
    double rAtt_h = k_att*Brs/3.0 - k_det_Rf*Brh;  // Rhodoferax hard <-> suspended
    double rAtt_m = k_att*Brs/3.0 - k_det_Rf*Brm;  // Rhodoferax medium
    double rAtt_e = k_att*Brs/3.0 - k_det_e*Bre;   // Rhodoferax easy
    double mAtt   = k_att*Bms     - k_det_Mt*Bma;  // methanogen attached <-> suspended
    double sAtt   = k_att*Bss     - k_det_SRB*Bsa; // sulfate reducer attached <-> suspended

    /* --- CARRYING-CAPACITY (logistic) GROWTH LIMITER for the 8 attached biofilm pools ---
     *  Real biofilms self-limit their density. Numerically, this is what keeps biomass from
     *  overshooting the CA density cap and thrashing the push-pull spreader (the "[CA] Stuck
     *  in push-pull loop" crash). GROWTH of the attached pools is throttled by
     *      f_cap = 1 - (total attached biomass)/B_cap        (clamped to [0,1])
     *  so biomass ASYMPTOTES to the cap instead of shooting past it -> the CA almost never
     *  has to spread, and when it does the excess is tiny and converges in a few sweeps.
     *  NOTE: only the GROWTH (yield x uptake) is throttled. The substrate reactions (U/Fe/S
     *  reduction, subsR above) are left at the full biomass-proportional rate, so cells at
     *  capacity keep catalysing at full speed - only their further growth stops (stationary
     *  phase / maintenance). Decay and attachment are unchanged. Suspended (planktonic) pools
     *  are NOT biofilm and are not capped.
     *  B_cap MUST equal the XML <maximum_biomass_density> (kg/m3); update both together. */
    constexpr double B_cap = 1.0e-4;                          // == XML maximum_biomass_density
    double Bfilm_tot = Bgh + Bgm + Bge + Brh + Brm + Bre + Bma + Bsa;   // total attached biofilm here
    double f_cap = 1.0 - Bfilm_tot / B_cap;
    if (f_cap < 0.0) f_cap = 0.0;                             // at/above the cap: no further growth (never negative)

    // Geobacter attached pools: growth throttled by f_cap; decay & attachment unchanged
    bioR[GEO_H] = f_cap*(Y_Geo_Fe*rGh_Fe + Y_Geo_U*rGh_U) - kd*Bgh       + gAtt_h;
    bioR[GEO_M] = f_cap*(Y_Geo_Fe*rGm_Fe + Y_Geo_U*rGm_U) - kd*Bgm       + gAtt_m;
    bioR[GEO_E] = f_cap*(Y_Geo_Fe*rGe_Fe + Y_Geo_U*rGe_U) - kd_Geo_e*Bge + gAtt_e;
    bioR[GEO_S] = Y_Geo_U*rGs_U   - kd*Bgs         - (gAtt_h+gAtt_m+gAtt_e);   // suspended: U only; NOT capped
    // Rhodoferax pools: grow on Fe(III) only (attached throttled by f_cap)
    bioR[RF_H]  = f_cap*(Y_Rf*rRh_Fe) - kd*Brh       + rAtt_h;
    bioR[RF_M]  = f_cap*(Y_Rf*rRm_Fe) - kd*Brm       + rAtt_m;
    bioR[RF_E]  = f_cap*(Y_Rf*rRe_Fe) - kd*Bre       + rAtt_e;
    bioR[RF_S]  =             - kd*Brs        - (rAtt_h+rAtt_m+rAtt_e);   // suspended: NOT capped
    // methanogen pools (attached throttled)
    bioR[MT_A]  = f_cap*(Y_Mt*rMa) - kd*Bma + mAtt;
    bioR[MT_S]  = Y_Mt*rMs - kd*Bms - mAtt;                              // suspended: NOT capped
    // sulfate reducer pools (attached throttled)
    bioR[SRB_A] = f_cap*(Y_SRB*rSa) - kd*Bsa + sAtt;
    bioR[SRB_S] = Y_SRB*rSs - kd*Bss - sAtt;                             // suspended: NOT capped

    // report using total Geobacter biomass and its net growth (reporting only)
    double Geo_tot = Bgh+Bgm+Bge+Bgs;
    double dGeo    = bioR[GEO_H]+bioR[GEO_M]+bioR[GEO_E]+bioR[GEO_S];
    KineticsStats::accumulate(Geo_tot, Ac, dGeo);

    // final safety net: if any rate came out not-a-number or infinite, reset it to 0 (T9 safeguard)
    for (size_t i = 0; i < subsR.size(); ++i)
        if (std::isnan(subsR[i]) || std::isinf(subsR[i])) subsR[i] = 0.0;
    for (size_t i = 0; i < bioR.size(); ++i)
        if (std::isnan(bioR[i]) || std::isinf(bioR[i])) bioR[i] = 0.0;
}

#endif // DEFINE_KINETICS_HH   // closes the include guard opened at the top
