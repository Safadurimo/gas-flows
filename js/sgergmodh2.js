// sgergmodh2.js  –  SGERG-mod-H2 compression factor (DVGW G 685-6:2024 Anhang H)
// Based on SGERG-88 (ISO 12213-3) with 6 modifications per DVGW G 685-6 (A) August 2024.
//
// The following 6 changes to the FORTRAN reference code (DIN EN ISO 12213-3:2010-01)
// are applied exactly as described in Anhang H of DVGW G 685-6 (A) August 2024, p. 36:
//
//   Change 1 (Fortran SGERG1 line 10, comment "REL. MASS OUT OF RANGE"):
//     Lower relative density bound reduced from 0.55 to 0.06.
//
//   Change 2 (Fortran SGERG1 line 12, comment "CALOR. VALUE OUT OF RANGE"):
//     Lower calorific value limit reduced from 20.0 to 6.0 MJ/m³.
//     (Note: this range check was absent from the base JS translation; added here
//      with the modified limit to reflect the Fortran behaviour.)
//
//   Change 3 (Fortran SGERG1 line 13, comment "CONFLICTING INPUT"):
//     Check removed.
//     (Note: this check was not present in the base JS translation of SGERG-88;
//      it is therefore already absent and no code change is required here.)
//
//   Change 4 (Fortran SGERG1 line 15):
//     X7 = X5*0.0964D0  →  X7 = 0D0
//     CO mole fraction is no longer derived from H2 content (Stadtgas legacy removed).
//
//   Change 5 (Fortran SGERG1 line 19):
//     BEFF = -0.065D0  →  BEFF = 0D0
//     Initial effective second virial coefficient no longer contains CO correction.
//
//   Change 6 (Fortran SGERG1 lines 52/53, comment "CONFLICTING RESULT FOR N2"):
//     Check removed.
//
// References:
//   DVGW G 685-6 (A) August 2024 – Abrechnungsbrennwert und Abrechnungszustandszahl,
//     Anhang H, p. 36.
//   DVGW Research Report DBI GUT 0083H:
//     https://www.dvgw.de/medien/dvgw/gas/infrastruktur/dvgw-pk-1-5-3-forschungsbericht-sgerg88-mod-h2-eng.pdf
//
// Variable names, subroutine names and structure kept identical to the FORTRAN
// original (and to sgerg.js) so that an auditor can verify changes by diff.
//
// COMMON-block state is passed as an object 's' between subroutines;
// all other names are unchanged.
//
// Units (same as FORTRAN):
//   HS  calorific value   MJ/m³  (metering at 0 °C, 1.01325 bar;
//                                  combustion reference 25 °C)
//   RM  relative density  –      (gas / air, both at 0 °C, 1.01325 bar)
//   P   pressure          bar
//   TC  temperature       °C
//   Z   compression factor –
//   D   molar density     mol/L  (= kmol/m³)

/* ─── BLOCK DATA ────────────────────────────────────────────────────────────
   Coefficients from the ISO 12213-3 / GERG-88 reference implementation.
   Unchanged from SGERG-88.                                                  */

// /BBLOK/
const BR11H0_M = [-0.425468,    0.286500e-2,  -0.462073e-5 ];
const BR11H1_M = [ 0.877118e-3,-0.556281e-5,   0.881510e-8 ];
const BR11H2_M = [-0.824747e-6, 0.431436e-8,  -0.608319e-11];
const BR22_M   = [-0.144600,    0.740910e-3,  -0.911950e-6 ];
const BR23_M   = [-0.339693,    0.161176e-2,  -0.204429e-5 ];
const BR33_M   = [-0.868340,    0.403760e-2,  -0.516570e-5 ];
const BR15_M   = [-0.521280e-1, 0.271570e-3,  -0.25e-6     ];
const BR17_M   = [-0.687290e-1,-0.239381e-5,   0.518195e-6 ];
const BR55_M   = [-0.110596e-2, 0.813385e-4,  -0.987220e-7 ];
const BR77_M   = [-0.130820,    0.602540e-3,  -0.644300e-6 ];
const B25_M    =  0.012;

// /CBLOK/
const CR111H0_M = [-0.302488,    0.195861e-2,  -0.316302e-5 ];
const CR111H1_M = [ 0.646422e-3,-0.422876e-5,   0.688157e-8 ];
const CR111H2_M = [-0.332805e-6, 0.223160e-8,  -0.367713e-11];
const CR222_M   = [ 0.784980e-2,-0.398950e-4,   0.611870e-7 ];
const CR223_M   = [ 0.552066e-2,-0.168609e-4,   0.157169e-7 ];
const CR233_M   = [ 0.358783e-2, 0.806674e-5,  -0.325798e-7 ];
const CR333_M   = [ 0.205130e-2, 0.348880e-4,  -0.837030e-7 ];
const CR555_M   = [ 0.104711e-2,-0.364887e-5,   0.467095e-8 ];
const CR117_M   = [ 0.736748e-2,-0.276578e-4,   0.343051e-7 ];

// /ZETA/
const Z12_M  =  0.72;
const Z13_M  = -0.865;
const Y12_M  =  0.92;
const Y13_M  =  0.92;
const Y123_M =  1.10;
const Y115_M =  1.2;

// /MBLOK/
const GM1R0_M = -2.709328;
const GM1R1_M =  0.021062199;
const GM2_M   = 28.0135;
const GM3_M   = 44.010;
const GM5_M   =  2.0159;
const GM7_M   = 28.010;
const FA_M    = 22.414097;
const T0_M    = 273.15;
const RL_M    =  1.292923;
const H5_M    = 285.83;
const H7_M    = 282.98;
const R_M     =  0.0831451;

/* ─── SUBROUTINE SMBER ───────────────────────────────────────────────────────
   Unchanged from SGERG-88.                                                   */
function SMBER_M(H, s) {
    const GM1 = GM1R0_M + GM1R1_M * H;
    s.X1 = (s.HS - (s.X5 * H5_M + s.X7 * H7_M) * s.AMOL) / (H * s.AMOL);
    s.X2 = 1.0 - s.X1 - s.X3 - s.X5 - s.X7;
    return (s.X1 * GM1 + s.X2 * GM2_M + s.X3 * GM3_M + s.X5 * GM5_M + s.X7 * GM7_M) * s.AMOL;
}

/* ─── SUBROUTINE B11BER ──────────────────────────────────────────────────────
   Unchanged from SGERG-88.                                                   */
function B11BER_M(T, H) {
    const T2 = T * T;
    return (BR11H0_M[0] + BR11H0_M[1]*T + BR11H0_M[2]*T2)
         + (BR11H1_M[0] + BR11H1_M[1]*T + BR11H1_M[2]*T2) * H
         + (BR11H2_M[0] + BR11H2_M[1]*T + BR11H2_M[2]*T2) * H * H;
}

/* ─── SUBROUTINE BBER ────────────────────────────────────────────────────────
   Unchanged from SGERG-88.                                                   */
function BBER_M(T, B11, s) {
    const T2  = T * T;
    const B22 = BR22_M[0] + BR22_M[1]*T + BR22_M[2]*T2;
    const B23 = BR23_M[0] + BR23_M[1]*T + BR23_M[2]*T2;
    const B33 = BR33_M[0] + BR33_M[1]*T + BR33_M[2]*T2;
    const B15 = BR15_M[0] + BR15_M[1]*T + BR15_M[2]*T2;
    const B55 = BR55_M[0] + BR55_M[1]*T + BR55_M[2]*T2;
    const B17 = BR17_M[0] + BR17_M[1]*T + BR17_M[2]*T2;
    const B77 = BR77_M[0] + BR77_M[1]*T + BR77_M[2]*T2;

    const BA13 = B11 * B33;
    if (BA13 < 0.0) throw new Error('BBER: no solution (BA13 < 0)');

    const ZZZ  = Z12_M + Math.pow(320.0 - T, 2) * 1.875e-5;

    return s.X11*B11  + s.X12*ZZZ*(B11+B22) + 2.0*s.X13*Z13_M*Math.sqrt(BA13)
         + s.X22*B22  + 2.0*s.X23*B23       + s.X33*B33   + s.X55*B55
         + 2.0*s.X15*B15 + 2.0*s.X25*B25_M  + 2.0*s.X17*B17 + s.X77*B77;
}

/* ─── SUBROUTINE CBER ────────────────────────────────────────────────────────
   Unchanged from SGERG-88.                                                   */
function CBER_M(T, H, s) {
    const T2 = T * T;
    const C111 = (CR111H0_M[0] + CR111H0_M[1]*T + CR111H0_M[2]*T2)
               + (CR111H1_M[0] + CR111H1_M[1]*T + CR111H1_M[2]*T2) * H
               + (CR111H2_M[0] + CR111H2_M[1]*T + CR111H2_M[2]*T2) * H * H;
    const C222 = CR222_M[0] + CR222_M[1]*T + CR222_M[2]*T2;
    const C223 = CR223_M[0] + CR223_M[1]*T + CR223_M[2]*T2;
    const C233 = CR233_M[0] + CR233_M[1]*T + CR233_M[2]*T2;
    const C333 = CR333_M[0] + CR333_M[1]*T + CR333_M[2]*T2;
    const C555 = CR555_M[0] + CR555_M[1]*T + CR555_M[2]*T2;
    const C117 = CR117_M[0] + CR117_M[1]*T + CR117_M[2]*T2;

    const CA112 = C111 * C111 * C222;
    const CA113 = C111 * C111 * C333;
    const CA122 = C111 * C222 * C222;
    const CA123 = C111 * C222 * C333;
    const CA133 = C111 * C333 * C333;
    const CA115 = C111 * C111 * C555;

    if (CA112 < 0.0 || CA113 < 0.0 || CA122 < 0.0
     || CA123 < 0.0 || CA133 < 0.0 || CA115 < 0.0)
        throw new Error('CBER: no solution (negative C-product)');

    const D3REP = 1.0 / 3.0;

    return  s.X1*s.X11       * C111
          + 3.0*s.X11*s.X2   * Math.pow(CA112, D3REP) * (Y12_M + (T - 270.0) * 0.0013)
          + 3.0*s.X11*s.X3   * Math.pow(CA113, D3REP) * Y13_M
          + 3.0*s.X1*s.X15   * Math.pow(CA115, D3REP) * Y115_M
          + 3.0*s.X1*s.X22   * Math.pow(CA122, D3REP) * (Y12_M + (T - 270.0) * 0.0013)
          + 6.0*s.X1*s.X2*s.X3 * Math.pow(CA123, D3REP) * Y123_M
          + 3.0*s.X1*s.X33   * Math.pow(CA133, D3REP) * Y13_M
          + s.X22*s.X2        * C222
          + 3.0*s.X22*s.X3   * C223
          + 3.0*s.X2*s.X33   * C233
          + s.X3*s.X33        * C333
          + s.X5*s.X55        * C555
          + 3.0*s.X11*s.X7   * C117;
}

/* ─── SUBROUTINE ITER ────────────────────────────────────────────────────────
   Unchanged from SGERG-88.                                                   */
function ITER_M(P, T, B, C) {
    const RT  = R_M * T;
    const RTP = RT / P;
    let V = RTP + B;
    for (let KK = 0; KK < 20; KK++) {
        V = RTP * (1.0 + B/V + C/(V*V));
        const Z  = 1.0 + B/V + C/(V*V);
        const PA = RT / V * Z;
        if (Math.abs(PA - P) < 1e-5) return { V, Z };
    }
    throw new Error('ITER: no convergence (#3)');
}

/* ─── SUBROUTINE SGERG1 ──────────────────────────────────────────────────────
   Modified per DVGW G 685-6 (A) August 2024, Anhang H.
   Changes vs SGERG-88 are annotated with "Anhang H change N".               */
function SGERG1_M(P, TC, Q3, Q5, QM, RM) {
    const s = {
        HS:   QM,
        X3:   Q3,
        X5:   Q5,
        X7:   0,            // Anhang H change 4 (Fortran line 15): was Q5 * 0.0964
        X1: 0, X2: 0,
        X11:0, X12:0, X13:0, X22:0, X23:0, X33: Q3*Q3,
        X55: Q5*Q5,  X77: 0, X15:0, X17:0, X25:0,
        AMOL: 0
    };
    s.X77 = s.X7 * s.X7;

    const SM = RM * RL_M;

    // Anhang H change 1 (Fortran line 10, "REL. MASS OUT OF RANGE"):
    // lower density bound reduced from 0.55 to 0.06
    if ((0.06 + 0.97*s.X3 - 0.45*s.X5) > RM)
        throw new Error('Relative mass out of range');
    // Anhang H change 3 (Fortran line 13, "CONFLICTING INPUT"):
    // check removed — was not present in the base JS translation

    // Anhang H change 5 (Fortran line 19): BEFF = -0.065D0  →  BEFF = 0D0
    let BEFF = 0;
    let H    = 1000.0;
    s.AMOL   = 1.0 / (FA_M + BEFF);

    let K  = 0;
    let KK = 0;

    outer: for (;;) {
        for (;;) {
            const SMT1 = SMBER_M(H, s);
            if (Math.abs(SM - SMT1) <= 1e-6) break;
            const SMT2 = SMBER_M(H + 1.0, s);
            H += (SM - SMT1) / (SMT2 - SMT1);
            KK++;
            if (KK > 20) throw new Error('No convergence #1');
        }

        s.X11 = s.X1 * s.X1;
        s.X12 = s.X1 * s.X2;
        s.X13 = s.X1 * s.X3;
        s.X22 = s.X2 * s.X2;
        s.X23 = s.X2 * s.X3;
        s.X25 = s.X2 * s.X5;
        s.X15 = s.X1 * s.X5;
        s.X17 = s.X1 * s.X7;

        const B11ref = B11BER_M(T0_M, H);
        BEFF   = BBER_M(T0_M, B11ref, s);
        s.AMOL = 1.0 / (FA_M + BEFF);

        const HSBER = s.X1 * H * s.AMOL + (s.X5 * H5_M + s.X7 * H7_M) * s.AMOL;
        if (Math.abs(s.HS - HSBER) <= 1e-4) break outer;
        K++;
        if (K > 20) throw new Error('No convergence #2');
    }

    // Validate calculated nitrogen fraction X2
    if (s.X2 < -0.01 || s.X2 > 0.5) throw new Error('Calculated N2 out of range');
    if (s.X2 + s.X3 > 0.5)          throw new Error('N2 + CO2 out of range');
    // Anhang H change 6 (Fortran lines 52/53, "CONFLICTING RESULT FOR N2"): check removed

    const T    = TC + T0_M;
    const B11  = B11BER_M(T, H);
    const B    = BBER_M(T, B11, s);
    const C    = CBER_M(T, H, s);
    const { V, Z } = ITER_M(P, T, B, C);
    const D    = 1.0 / V;

    return { Z, D, X2: s.X2 };
}

/* ─── PUBLIC ENTRY POINT ─────────────────────────────────────────────────────
   Mirrors SUBROUTINE SGERG(X2, X3, HS, RM, X5, P, TC, Z, D)

   Inputs
     X3   CO2 mole fraction          (0.0 – 0.20, extended to 0.30)
     HS   calorific value  [MJ/m³]   (6   – 48,  metering 0 °C / 1.01325 bar)
     RM   relative density           (0.06 – 0.90, relative to air)
     X5   H2  mole fraction          (0.0 – 0.30)
     P    pressure  [bar]            (0   – 120)
     TC   temperature  [°C]          (-23 – 65)

   Returns { Z, D [mol/L = kmol/m³], X2 [calculated N2 mole fraction] }

   Application ranges per DVGW G 685-6 (A) August 2024, Tabelle B.4.          */
function sgergmodh2(X3, HS, RM, X5, P, TC) {
    if (P  <  0.0 || P  > 120.0) throw new Error('Pressure out of range (0 – 120 bar)');
    if (TC < -23.0 || TC >  65.0) throw new Error('Temperature out of range (−23 – 65 °C)');
    // Anhang H change 2 (Fortran line 12, "CALOR. VALUE OUT OF RANGE"):
    // lower calorific value limit reduced from 20.0 to 6.0 MJ/m³
    if (HS < 6.0 || HS > 48.0) throw new Error('Calorific value out of range (6 – 48 MJ/m³)');
    return SGERG1_M(P, TC, X3, X5, HS, RM);
}

window.sgergmodh2 = sgergmodh2;
