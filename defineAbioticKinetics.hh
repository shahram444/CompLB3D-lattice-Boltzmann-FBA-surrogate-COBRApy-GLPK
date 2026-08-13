/* ============================================================================
 * defineAbioticKinetics.hh  —  FULL ZHAO NETWORK : ABIOTIC (CHEMICAL) REACTIONS
 *                              R108, R112, R113, R114   (R111 = A1 is OFF)
 * ============================================================================
 *
 * HOW TO READ THIS FILE (for readers who do not program):
 *
 *   - This is the "recipe" file for the purely chemical reactions - the ones
 *     that happen without microbes. The solver hands it the local chemical
 *     amounts in one voxel and it hands back how fast each chemical is changing
 *     there. It is called once per voxel, every time step. Pair it with the
 *     biotic file defineKinetics.hh (reactions R74-R107).
 *
 *   - Anything after "//" or between the block-comment markers is a note for
 *     humans; the computer ignores it.
 *
 *   - The 14 reacting chemicals C0..C13 are read by position, counting from 0:
 *        C0=acetate  C1=U(VI)   C2=U(IV)    C3=Fe(II)
 *        C4=Fe(III) easy   C5=Fe(III) medium   C6=Fe(III) hard
 *        C7=sulfate  C8=sulfide  C9=FeS(mackinawite)
 *        C10=bicarbonate  C11=H+ (proton)  C12=methane  C13=S0 (elemental sulfur, made here)
 *     (The solver's chemical vector also carries 10 equilibrium-only components
 *      C14..C23 - major ions and mineral surface sites - which this file never touches.)
 *
 *   THE FOUR CHEMICAL REACTIONS HERE (word form; solids marked (s)):
 *     R108  U(VI) reduced by sulfide:    U(VI) + sulfide -> U(IV)(s) + S0(s) + H+
 *     R112  FeS precipitation (A2):      Fe(II) + sulfide -> FeS(s) + H+
 *     R113  goethite reduced by sulfide: 2 Fe(III)hard + sulfide + 5 H+ -> 2 Fe(II) + S0(s)
 *     R114  U(VI) reduced by FeS (A3):   FeS(s) + 4 U(VI) + 8 H+ -> 4 U(IV)(s) + Fe(II) + sulfate
 *
 *   R111 (A1: dissolved Fe(II) + U(VI) -> Fe(III) + U(IV)) is deliberately OFF:
 *     Zhao's thermodynamic analysis (paper Figs. 4-5) shows this reaction has
 *     Gibbs energy > 0 below ~2 mM Fe(II)/U(VI), so it does not proceed in the
 *     field. It is documented here but never computed.
 *
 *   R109-R110 (sorbed U(VI) reduced by sorbed Fe(II) on iron-oxide / other surfaces)
 *     are NOT coded as separate terms. Following Zhao and the Option-A guide, their
 *     effect is carried by the same FeS surface-fraction factor frac_ab used in R114:
 *     as FeS covers the surface, U(VI) reduction shifts onto the mineral pathway.
 *     So all of Zhao's kinetic reactions are accounted for: R74-R107 (biotic file),
 *     R108 and R112-R114 here, R111 off, and R109-R110 via the surface fraction.
 *
 *   THE KEY COMPETITION (Zhao Eq. 2). R114 is throttled by frac_ab, the share of
 *   the surface now COVERED by FeS: frac_ab = FeS/(FeS + total Fe(III) + eps).
 *   It is the exact complement of the biotic file's frac_bio, so as FeS builds up
 *   the U(VI) reduction hands off from the biological pathway to this chemical one.
 *
 *   ELEMENT BOOKKEEPING (why the balancing terms below). R114 consumes one FeS per
 *   four U(VI) and releases the FeS's iron as Fe(II) and its sulfur as sulfate
 *   (Zhao's full stoichiometry). Writing those release terms is what makes iron
 *   and sulfur close exactly - it is the fix for the small iron/sulfur imbalance
 *   that the simplified T6/T9 set had (which dropped the FeS-oxidation products).
 *
 *   SAFEGUARDS CARRIED OVER FROM TEST T9: every amount floored at zero before use;
 *   a per-step stability limit scales reactions so no chemical is over-drawn in a
 *   step; a final not-a-number/infinity guard. Rate constants are PER-SECOND (Zhao per-day / 86400), so the
 *   solver advances them with ade_dt in SECONDS (see the biotic file's
 *   note on dt_kinetics and the full-network docs, Section 12).
 * ============================================================================ */
#ifndef DEFINE_ABIOTIC_KINETICS_HH   // include guard: read this file only once at build time
#define DEFINE_ABIOTIC_KINETICS_HH
/* ===========================================================================
 * [TRUE-UNITS 2026-07-24] ABIOTIC RATE CONSTANTS CONVERTED TO PHYSICAL PER-SECOND
 * (per-DAY / 86400): k_FeS,k_A4 in L.mol^-1.s^-1 ; K_U_FeS in s^-1 ;
 * K_U_HS in (mmol/L)^-0.54.s^-1 (n_HS=0.54 dimensionless, unchanged).
 * dt_kinetics set to the run ade_dt (~3.0e-4 s).
 * ============================================================================ */


#include <vector>       // the "list of numbers" type used to pass amounts in and rates out
#include <cmath>        // math helpers: std::pow (a^b) for R108, and the not-a-number test
#include <algorithm>    // std::max ("take the larger of two numbers")
#include <iostream>     // screen printing, available if a message is ever needed

/* ---------------------------------------------------------------------------
 * CHEMICAL POSITIONS (names for the list positions above)
 * --------------------------------------------------------------------------- */
namespace ZhaoAbioIdx {
    constexpr int AC=0, U6=1, U4=2, FE2=3, FE3E=4, FE3M=5, FE3H=6,
                  SO4=7, HS=8, FES=9, HCO3=10, HP=11, CH4=12, S0=13;
}

/* ---------------------------------------------------------------------------
 * THE FIXED NUMBERS (reaction constants) - Zhao Supplementary Table S3.
 *   Rate constants are PER-SECOND (Zhao per-day / 86400). "constexpr double NAME = VALUE;"
 *   gives a fixed decimal number a NAME so we can use it below.
 * --------------------------------------------------------------------------- */
namespace AbioticParams {
    constexpr double k_FeS   = 1.5;  // FeS precipitation rate constant, R112 (L/mol/s; Rickard 1995).
                                          // NOTE: Zhao's SI argues the subsurface value is ~1/10 of this
                                          // mixing-chamber number; set to 0.15 for the field-rate case.
    constexpr double K_U_FeS = 4.63e-6;       // U(VI) reduction by FeS, R114 (per second; Veeramani 2013 ES&T 47:2361, Zhao Eq.2)
    constexpr double k_A4    = 2.41e-4;      // goethite (hard Fe) reduction by sulfide, R113 (L/mol/s; Poulton 2004 / Fang 2009)
    constexpr double K_U_HS  = 7.16e-3;     // U(VI) reduction by sulfide, R108 (Hua 2006). = 618.3/day / 86400 (pure time conversion, preserves original model behavior). OPEN QUESTION: sulfide concentration basis unverified - if Hua's pre-factor 0.0103 (mol/L)^-0.54 min^-1 applies with C[HS] in mol/L, this should be 1.72e-4 (~42x lower). Verify vs Hua 2006 full text + model SI with Dr. Meile before adjusting.
    constexpr double n_HS    = 0.54;      // fractional reaction order in sulfide for R108 (Hua 2006)

    constexpr double eps      = 1.0e-12;  // tiny number so we never divide by zero
    constexpr double MAX_RATE_FRACTION = 0.5; // stability limit: a reaction may use at most half of a chemical per step
    // TIME UNITS: the solver advances these rates with ade_dt IN SECONDS = ((tau-0.5)/3)*dx^2/D0
    // (D0 = acetate pore diffusion 1e-9). dx=1 um -> 1.0e-4 s ; dx=2 um -> 4.0e-4 s. dt_kinetics is used
    // ONLY by the stability limit and MUST equal the run's ade_dt (see "[ADE] dt" in the log).
    constexpr double dt_kinetics = 1.0e-4;    // [seconds] = the run's printed [ADE] dt for this dx=1um duct (~1.0e-4). Set equal to ade_dt for any geometry.

    // -- Zhao stoichiometric coefficients (Supplementary Table S1) --
    constexpr double s_Fe3_per_HS = 2.0;  // mol hard-Fe(III) reduced per mol sulfide in R113
    constexpr double s_Hp_R113    = 5.0;  // protons consumed per mol sulfide in R113
    constexpr double s_U_per_FeS  = 4.0;  // mol U(VI) reduced per mol FeS in R114
    constexpr double s_Hp_R114    = 8.0;  // protons consumed per 4 U(VI) (= 2 per U(VI)) in R114
}

/* ---------------------------------------------------------------------------
 * SIMPLE RUN STATISTICS  (reporting only; no effect on the chemistry)
 * --------------------------------------------------------------------------- */
namespace AbioticKineticsStats {
    static double iter_total_reaction = 0.0;   // running total of reaction activity this step
    static long   iter_cells_reacting = 0;     // running count of cells that reacted this step
    inline void resetIteration() { iter_total_reaction = 0.0; iter_cells_reacting = 0; }
    inline void accumulate(double rate) {
        if (std::abs(rate) > 1e-20) { iter_cells_reacting++; iter_total_reaction += rate; }
    }
}

/* ---------------------------------------------------------------------------
 * THE MAIN ABIOTIC REACTION ROUTINE
 *   Runs once per voxel per step. IN: C (the chemical amounts). OUT: subsR (how
 *   fast each chemical changes). The "&" means the solver reads subsR back after
 *   the routine finishes. mask labels the voxel type (not used here). This file
 *   reads/writes only the 14 reacting chemicals C0..C13.
 * --------------------------------------------------------------------------- */
void defineAbioticRxnKinetics(
    std::vector<double> C,          // IN : the chemical amounts in this voxel (reacting species C0..C13; see guide)
    std::vector<double>& subsR,     // OUT: how fast each chemical changes (this file writes C0..C13)
    plb::plint mask                 // IN : voxel label; not used here (the processor decides where to call)
) {
    using namespace AbioticParams;
    using namespace ZhaoAbioIdx;

    for (size_t i = 0; i < subsR.size(); ++i) subsR[i] = 0.0;   // start every chemical rate at 0
    if (C.size() < 14) return;                                  // if the list is short, stop and do nothing

    // read the chemicals; std::max(x,0.0) floors any tiny negative at zero (T9 safeguard)
    double cU6  = std::max(C[U6],   0.0);   // dissolved U(VI)   (local name differs from the index U6)
    double Fe2  = std::max(C[FE2],  0.0);   // dissolved Fe(II)
    double Fe3h = std::max(C[FE3H], 0.0);   // solid Fe(III), hard class (goethite; the R113 reactant)
    double cHS  = std::max(C[HS],   0.0);   // dissolved sulfide (free HS- bucket for now; totalled below)
    double FeS  = std::max(C[FES],  0.0);   // FeS mineral present so far

    // === FULL-NETWORK TOTALS (kinetics_on_totals; NO change needed anywhere else) ======
    //  Same idea as the biotic file: with the full equilibrium ON, C[U6]/C[FE2]/C[HS] hold
    //  only the FREE ion. The abiotic rates (R108,R112,R114) are Zhao's, written on TOTAL
    //  dissolved concentrations, so rebuild them here from the AQUEOUS complexes (immobile
    //  surface-bound complexes excluded). Indices fixed by the 95-substrate layout. Guarded
    //  by C.size()>24 so the SAME file also runs unchanged with equilibrium OFF (24 subs).
    if (C.size() > 24) {
        cU6 += 2*(C[65]+C[66]+C[67]) + 3*(C[68]+C[69]+C[70]+C[71])
             + C[72]+C[73]+C[74]+C[75]+C[76]+C[77]+C[78]+C[79]+C[80];   // aqueous U(VI) complexes
        Fe2 += C[34]+C[35]+C[36]+C[37]+C[38]+C[39]+C[40]+C[41]+C[83];   // aqueous Fe(II) complexes (+ FeAc)
        cHS += C[43]+C[64];                                            // H2S + S2-
    }
    // ================================================================================
    double Fe3tot = std::max(C[FE3E],0.0) + std::max(C[FE3M],0.0) + Fe3h;  // total Fe(III) surface

    // -- reaction rates (per second) --
    double R108 = K_U_HS * cU6 * std::pow(cHS, n_HS); // R108: U(VI) reduced by sulfide (fractional order n in HS)
    double R112 = k_FeS  * Fe2 * cHS;                 // R112: FeS forms from Fe(II) and sulfide
    double R113 = k_A4   * Fe3h * cHS;                // R113: sulfide reduces hard Fe(III) (per mol sulfide)
    double frac_ab = FeS / (FeS + Fe3tot + eps);      // share of surface covered by FeS (0 -> 1 as FeS builds)
    double R114 = K_U_FeS * frac_ab * cU6;            // R114: FeS-surface reduction of U(VI) (Zhao Eq.2; U6-consumption rate)

    /* --- PER-STEP STABILITY LIMIT (T9 safeguard) ---
     *   Trim each rate so no chemical loses more than half of what it holds in a
     *   step. Sulfide is shared by R108, R112, R113, so if together they want too
     *   much sulfide, shrink those three by one common factor (keeps balances exact). */
    double capFe2 = Fe2  * MAX_RATE_FRACTION / dt_kinetics;         // most Fe(II) R112 may use
    if (R112 > capFe2) R112 = capFe2;
    double capFeh = Fe3h * MAX_RATE_FRACTION / dt_kinetics;         // most hard Fe(III) R113 may use (2 per sulfide)
    if (s_Fe3_per_HS*R113 > capFeh && R113 > 0.0) R113 = capFeh / s_Fe3_per_HS;
    double hsDraw = R108 + R112 + R113;                            // total sulfide demanded by the three
    double capHS  = cHS  * MAX_RATE_FRACTION / dt_kinetics;
    if (hsDraw > capHS && hsDraw > 0.0) { double s = capHS/hsDraw; R108*=s; R112*=s; R113*=s; }
    double uDraw  = R108 + R114;                                   // total U(VI) demanded by R108 + R114
    double capU   = cU6  * MAX_RATE_FRACTION / dt_kinetics;
    if (uDraw > capU && uDraw > 0.0) { double s = capU/uDraw; R108*=s; R114*=s; }
    double capFeS = FeS  * MAX_RATE_FRACTION / dt_kinetics;         // most FeS R114 may consume (1 per 4 U(VI))
    if (R114/s_U_per_FeS > capFeS && R114 > 0.0) R114 = capFeS * s_U_per_FeS;

    /* --- WRITE HOW FAST EACH CHEMICAL CHANGES (with the balancing terms that close U, Fe, S) --- */
    // helper: remove a total consumption 'cons' of a dissolved element spread over its buckets
    // in proportion to content (conserves the element, no bucket driven negative even when the
    // free ion is a tiny fraction of the total). Collapses to the single free bucket when eq is OFF.
    auto drawFrom = [&](double cons, double tot, std::initializer_list<int> idx){
        if (tot <= 0.0 || cons <= 0.0) return;
        for (int i : idx) if (i < (int)C.size()) subsR[i] += -cons * std::max(C[i],0.0) / tot;
    };

    // uranium: consumed by R108 and R114 (spread over free + aqueous U complexes); same amount -> immobile U(IV)
    drawFrom(R108 + R114, cU6, {U6,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80});
    subsR[U4] = +R108 + R114;                          // U(IV) up (uranium conserved)
    // iron: R113 releases 2 Fe(II) per sulfide from hard Fe(III); R112 locks Fe(II) into FeS;
    //       R114 releases the FeS's iron back as Fe(II) (1 Fe per 4 U(VI))
    subsR[FE3H] = -s_Fe3_per_HS * R113;                // hard Fe(III) consumed by R113
    drawFrom(R112, Fe2, {FE2,34,35,36,37,38,39,40,41,83});   // Fe(II) consumed by R112 (spread over free + Fe complexes)
    subsR[FE2] += s_Fe3_per_HS*R113 + R114/s_U_per_FeS;      // Fe(II) released by R113 and by FeS dissolution (to free pool)
    subsR[FES]  = +R112 - R114/s_U_per_FeS;            // FeS: made by R112, consumed by R114 (1 per 4 U(VI))
    // sulfur: R112 and R113 consume sulfide; R114 releases the FeS's sulfur as sulfate (1 per 4 U(VI));
    //         R108 and R113 turn some sulfide into elemental sulfur S0(s), tracked in C13 so sulfur closes
    drawFrom(R108 + R112 + R113, cHS, {HS,43,64});     // sulfide consumed (spread over free + H2S + S2-)
    subsR[SO4] = +R114/s_U_per_FeS;                    // sulfate released as FeS sulfur is oxidized in R114
    subsR[S0]  = +R108 + R113;                         // elemental sulfur produced by R108 and R113 (immobile)
    // protons: R108 and R112 release H+, R113 and R114 consume H+
    subsR[HP]  = +R108 + R112 - s_Hp_R113*R113 - (s_Hp_R114/s_U_per_FeS)*R114;

    AbioticKineticsStats::accumulate(R108 + R112 + R113 + R114);   // progress tally (reporting only)

    // final safety net: reset any not-a-number or infinite value to 0 (T9 safeguard)
    for (size_t i = 0; i < subsR.size(); ++i)
        if (std::isnan(subsR[i]) || std::isinf(subsR[i])) subsR[i] = 0.0;
}


/* ================================================================================================
 * ==============================  MINERAL DISSOLUTION RATE  ======================================
 * ================================================================================================
 *
 *  This is the dissolution counterpart of the precipitation reactions above, and like them it is
 *  YOUR chemistry: dissolutionVOP.hh only handles the geometry (which voxels are wetted, where the
 *  products go, when a voxel reopens).  It calls this function once per wetted mineral voxel per
 *  time step.
 *
 *  You only need to touch this if <dissolution><enabled>true</enabled></dissolution> is set in
 *  CompLaB.xml.  With dissolution off it is never called, and the empty default below is harmless.
 *
 *  ------------------------------------------------------------------------------------------------
 *  WHAT YOU ARE GIVEN
 *
 *    phaseId   which solid phase this voxel is, numbered 1, 2, 3 ... in the order the <phaseN>
 *              blocks appear in CompLaB.xml.  Branch on it to give each mineral its own rate law.
 *
 *    C         the water IN CONTACT with the mineral: the solute concentrations of the open
 *              neighbouring voxels, averaged.  Same substrate order as everywhere else, mol/L.
 *              This is not the mineral voxel's own water -- a solid voxel has none that moves.
 *
 *    mineral   how much of this mineral is left in the voxel, mol/L.  A voxel that is completely
 *              full holds <full_density>; one that is nearly gone holds almost nothing.  Use it
 *              if your rate depends on how much surface is left.
 *
 *    mask      the voxel's material number, if you want the rate to differ by region.
 *
 *  WHAT YOU MUST RETURN
 *
 *    mineralR  how fast the mineral disappears, mol/L/s.  MUST BE NEGATIVE (or zero for no
 *              reaction) -- this function only dissolves.  A positive value is ignored, because
 *              growing the mineral is what defineAbioticRxnKinetics() above is for.
 *
 *    subsR     what that releases into the water, mol/L/s, one per substrate.  POSITIVE means
 *              produced.  Do not write the mineral's own entry; the solver handles it from
 *              mineralR.
 *
 *  UNITS: per second, like every other rate in this file.  A published rate constant given per day
 *  must be divided by 86400 before it goes here.  That single mistake is the most common cause of a
 *  run that does nothing or explodes.
 *
 *  STABILITY: the solver will never let one step remove more than half the mineral present, and it
 *  scales everything you return by the SAME factor when it does, so your element balance survives
 *  the limiter exactly.  You do not need to add your own cap.
 *
 *  ------------------------------------------------------------------------------------------------
 *  WORKED EXAMPLE (commented out -- the shipped default does nothing)
 *
 *  Mackinawite dissolving in acid,   FeS + H+  ->  Fe2+ + HS-
 *  A simple first-order form, proportional to the proton concentration and to how much mineral is
 *  left (a crude stand-in for reactive surface area):
 *
 *      constexpr double k_FeS_diss = 1.0e-6;      // L/mol/s, per mol/L of remaining mineral
 *      if (phaseId == 1) {
 *          const double rate = k_FeS_diss * C[HP] * mineral;   // mol/L/s, positive magnitude
 *          mineralR   = -rate;                                 // mineral disappears
 *          subsR[FE2] = +rate;                                 // one Fe2+ per FeS
 *          subsR[HS]  = +rate;                                 // one HS- per FeS
 *          subsR[HP]  = -rate;                                 // one proton consumed
 *      }
 *
 *  Note the stoichiometry closes by construction: one mole of FeS removed produces exactly one mole
 *  each of Fe2+ and HS-.  Keep that property and the mass balance takes care of itself.
 *
 *  A saturation-dependent form -- rate proportional to how far from equilibrium the water is --
 *  would be more physical, but the equilibrium solver does not currently compute solubility
 *  products, so it cannot supply the saturation state.  First-order in the reactant is the
 *  practical choice today.
 * ================================================================================================
 */
inline void defineDissolutionRate(plb::plint phaseId,
                                  const std::vector<double>& C,
                                  double mineral,
                                  std::vector<double>& subsR,
                                  double& mineralR,
                                  plb::plint mask)
{
    (void) phaseId; (void) C; (void) mineral; (void) mask;

    mineralR = 0.0;                       // no dissolution by default
    for (size_t i = 0; i < subsR.size(); ++i) subsR[i] = 0.0;

    // ------------------------------------------------------------------------------------------
    // WRITE YOUR DISSOLUTION CHEMISTRY HERE.  See the worked example above.
    // ------------------------------------------------------------------------------------------

    // final safety net, matching defineAbioticRxnKinetics above
    if (std::isnan(mineralR) || std::isinf(mineralR)) mineralR = 0.0;
    for (size_t i = 0; i < subsR.size(); ++i)
        if (std::isnan(subsR[i]) || std::isinf(subsR[i])) subsR[i] = 0.0;
}

#endif // DEFINE_ABIOTIC_KINETICS_HH   // closes the include guard opened at the top
