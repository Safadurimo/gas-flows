// sgerg.js  –  SGERG-88 compression factor (ISO 12213-3)
// Translated from FORTRAN reference code by M. Jaeschke & J. Sikora, 07.08.96
// Variable names, subroutine names and structure kept identical to the FORTRAN
// original so that an auditor familiar with the ISO source can follow the code.
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
   Coefficients from the ISO 12213-3 / GERG-88 reference implementation.    */

// /BBLOK/
const BR11H0 = [-0.425468,    0.286500e-2,  -0.462073e-5 ];
const BR11H1 = [ 0.877118e-3,-0.556281e-5,   0.881510e-8 ];
const BR11H2 = [-0.824747e-6, 0.431436e-8,  -0.608319e-11];
const BR22   = [-0.144600,    0.740910e-3,  -0.911950e-6 ];
const BR23   = [-0.339693,    0.161176e-2,  -0.204429e-5 ];
const BR33   = [-0.868340,    0.403760e-2,  -0.516570e-5 ];
const BR15   = [-0.521280e-1, 0.271570e-3,  -0.25e-6     ];
const BR17   = [-0.687290e-1,-0.239381e-5,   0.518195e-6 ];
const BR55   = [-0.110596e-2, 0.813385e-4,  -0.987220e-7 ];
const BR77   = [-0.130820,    0.602540e-3,  -0.644300e-6 ];
const B25    =  0.012;

// /CBLOK/
const CR111H0 = [-0.302488,    0.195861e-2,  -0.316302e-5 ];
const CR111H1 = [ 0.646422e-3,-0.422876e-5,   0.688157e-8 ];
const CR111H2 = [-0.332805e-6, 0.223160e-8,  -0.367713e-11];
const CR222   = [ 0.784980e-2,-0.398950e-4,   0.611870e-7 ];
const CR223   = [ 0.552066e-2,-0.168609e-4,   0.157169e-7 ];
const CR233   = [ 0.358783e-2, 0.806674e-5,  -0.325798e-7 ];
const CR333   = [ 0.205130e-2, 0.348880e-4,  -0.837030e-7 ];
const CR555   = [ 0.104711e-2,-0.364887e-5,   0.467095e-8 ];
const CR117   = [ 0.736748e-2,-0.276578e-4,   0.343051e-7 ];

// /ZETA/
const Z12  =  0.72;
const Z13  = -0.865;
const Y12  =  0.92;
const Y13  =  0.92;
const Y123 =  1.10;
const Y115 =  1.2;

// /MBLOK/
const GM1R0 = -2.709328;       // molar-mass coefficient for CH4-rich fraction
const GM1R1 =  0.021062199;
const GM2   = 28.0135;         // molar mass N2  [g/mol]
const GM3   = 44.010;          // molar mass CO2 [g/mol]
const GM5   =  2.0159;         // molar mass H2  [g/mol]
const GM7   = 28.010;          // molar mass CO  [g/mol]
const FA    = 22.414097;       // ideal molar volume at T0, P0  [L/mol]
const T0    = 273.15;          // reference temperature  [K]
const RL    =  1.292923;       // air density at T0, P0  [kg/m³]
const H5    = 285.83;          // calorific value H2  [kJ/mol]
const H7    = 282.98;          // calorific value CO  [kJ/mol]
const R     =  0.0831451;      // gas constant  [L·bar / (mol·K)]

/* ─── SUBROUTINE SMBER ───────────────────────────────────────────────────────
   Calculates SM (mass density at reference conditions) for a given H.
   Side-effects: updates s.X1, s.X2  (like Fortran COMMON /XBLOK/).         */
function SMBER(H, s) {
    const GM1 = GM1R0 + GM1R1 * H;
    s.X1 = (s.HS - (s.X5 * H5 + s.X7 * H7) * s.AMOL) / (H * s.AMOL);
    s.X2 = 1.0 - s.X1 - s.X3 - s.X5 - s.X7;
    return (s.X1 * GM1 + s.X2 * GM2 + s.X3 * GM3 + s.X5 * GM5 + s.X7 * GM7) * s.AMOL;
}

/* ─── SUBROUTINE B11BER ──────────────────────────────────────────────────────
   Second virial coefficient for the hydrocarbon fraction.                   */
function B11BER(T, H) {
    const T2 = T * T;
    return (BR11H0[0] + BR11H0[1]*T + BR11H0[2]*T2)
         + (BR11H1[0] + BR11H1[1]*T + BR11H1[2]*T2) * H
         + (BR11H2[0] + BR11H2[1]*T + BR11H2[2]*T2) * H * H;
}

/* ─── SUBROUTINE BBER ────────────────────────────────────────────────────────
   Effective second virial coefficient of the mixture.                       */
function BBER(T, B11, s) {
    const T2  = T * T;
    const B22 = BR22[0] + BR22[1]*T + BR22[2]*T2;
    const B23 = BR23[0] + BR23[1]*T + BR23[2]*T2;
    const B33 = BR33[0] + BR33[1]*T + BR33[2]*T2;
    const B15 = BR15[0] + BR15[1]*T + BR15[2]*T2;
    const B55 = BR55[0] + BR55[1]*T + BR55[2]*T2;
    const B17 = BR17[0] + BR17[1]*T + BR17[2]*T2;
    const B77 = BR77[0] + BR77[1]*T + BR77[2]*T2;

    const BA13 = B11 * B33;
    if (BA13 < 0.0) throw new Error('BBER: no solution (BA13 < 0)');

    const ZZZ  = Z12 + Math.pow(320.0 - T, 2) * 1.875e-5;

    return s.X11*B11  + s.X12*ZZZ*(B11+B22) + 2.0*s.X13*Z13*Math.sqrt(BA13)
         + s.X22*B22  + 2.0*s.X23*B23       + s.X33*B33   + s.X55*B55
         + 2.0*s.X15*B15 + 2.0*s.X25*B25    + 2.0*s.X17*B17 + s.X77*B77;
}

/* ─── SUBROUTINE CBER ────────────────────────────────────────────────────────
   Effective third virial coefficient of the mixture.                        */
function CBER(T, H, s) {
    const T2 = T * T;
    const C111 = (CR111H0[0] + CR111H0[1]*T + CR111H0[2]*T2)
               + (CR111H1[0] + CR111H1[1]*T + CR111H1[2]*T2) * H
               + (CR111H2[0] + CR111H2[1]*T + CR111H2[2]*T2) * H * H;
    const C222 = CR222[0] + CR222[1]*T + CR222[2]*T2;
    const C223 = CR223[0] + CR223[1]*T + CR223[2]*T2;
    const C233 = CR233[0] + CR233[1]*T + CR233[2]*T2;
    const C333 = CR333[0] + CR333[1]*T + CR333[2]*T2;
    const C555 = CR555[0] + CR555[1]*T + CR555[2]*T2;
    const C117 = CR117[0] + CR117[1]*T + CR117[2]*T2;

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
          + 3.0*s.X11*s.X2   * Math.pow(CA112, D3REP) * (Y12 + (T - 270.0) * 0.0013)
          + 3.0*s.X11*s.X3   * Math.pow(CA113, D3REP) * Y13
          + 3.0*s.X1*s.X15   * Math.pow(CA115, D3REP) * Y115
          + 3.0*s.X1*s.X22   * Math.pow(CA122, D3REP) * (Y12 + (T - 270.0) * 0.0013)
          + 6.0*s.X1*s.X2*s.X3 * Math.pow(CA123, D3REP) * Y123
          + 3.0*s.X1*s.X33   * Math.pow(CA133, D3REP) * Y13
          + s.X22*s.X2        * C222
          + 3.0*s.X22*s.X3   * C223
          + 3.0*s.X2*s.X33   * C233
          + s.X3*s.X33        * C333
          + s.X5*s.X55        * C555
          + 3.0*s.X11*s.X7   * C117;
}

/* ─── SUBROUTINE ITER ────────────────────────────────────────────────────────
   Iterative solver for molar volume V and compression factor Z.             */
function ITER(P, T, B, C) {
    const RT  = R * T;
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
   Main calculation routine (SGERG1 in FORTRAN).                             */
function SGERG1(P, TC, Q3, Q5, QM, RM) {
    // State object simulates Fortran COMMON blocks /RBLOK/ /XBLOK/
    const s = {
        HS:   QM,           // calorific value  [MJ/m³]
        X3:   Q3,           // CO2 mole fraction
        X5:   Q5,           // H2  mole fraction
        X7:   Q5 * 0.0964,  // CO  mole fraction (calculated from H2)
        X1: 0, X2: 0,
        X11:0, X12:0, X13:0, X22:0, X23:0, X33: Q3*Q3,
        X55: Q5*Q5,  X77: 0, X15:0, X17:0, X25:0,
        AMOL: 0
    };
    s.X77 = s.X7 * s.X7;

    const SM = RM * RL;   // target mass density at reference conditions [kg/m³]

    if ((0.55 + 0.97*s.X3 - 0.45*s.X5) > RM)
        throw new Error('Conflicting input parameters');

    let BEFF = -0.065;
    let H    = 1000.0;
    s.AMOL   = 1.0 / (FA + BEFF);

    let K  = 0;
    let KK = 0;

    // ── combined inner/outer iteration (label 1 in FORTRAN) ──────────────
    outer: for (;;) {
        // Inner loop: find H such that SMBER(H) = SM  (secant method)
        for (;;) {
            const SMT1 = SMBER(H, s);
            if (Math.abs(SM - SMT1) <= 1e-6) break;
            const SMT2 = SMBER(H + 1.0, s);
            H += (SM - SMT1) / (SMT2 - SMT1);
            KK++;
            if (KK > 20) throw new Error('No convergence #1');
        }

        // Update cross-products (X11, X12, … after X1, X2 are known)
        s.X11 = s.X1 * s.X1;
        s.X12 = s.X1 * s.X2;
        s.X13 = s.X1 * s.X3;
        s.X22 = s.X2 * s.X2;
        s.X23 = s.X2 * s.X3;
        s.X25 = s.X2 * s.X5;
        s.X15 = s.X1 * s.X5;
        s.X17 = s.X1 * s.X7;

        // Update BEFF and AMOL at reference temperature T0
        const B11ref = B11BER(T0, H);
        BEFF   = BBER(T0, B11ref, s);
        s.AMOL = 1.0 / (FA + BEFF);

        // Outer convergence: calorific value match
        const HSBER = s.X1 * H * s.AMOL + (s.X5 * H5 + s.X7 * H7) * s.AMOL;
        if (Math.abs(s.HS - HSBER) <= 1e-4) break outer;
        K++;
        if (K > 20) throw new Error('No convergence #2');
    }

    // Validate calculated nitrogen fraction X2
    if (s.X2 < -0.01 || s.X2 > 0.5) throw new Error('Calculated N2 out of range');
    if (s.X2 + s.X3 > 0.5)          throw new Error('N2 + CO2 out of range');
    if ((0.55 + 0.4*s.X2 + 0.97*s.X3 - 0.45*s.X5) > RM)
        throw new Error('Conflicting result for N2');

    // Final calculation at operating conditions
    const T    = TC + T0;
    const B11  = B11BER(T, H);
    const B    = BBER(T, B11, s);
    const C    = CBER(T, H, s);
    const { V, Z } = ITER(P, T, B, C);
    const D    = 1.0 / V;

    return { Z, D, X2: s.X2 };
}

/* ─── PUBLIC ENTRY POINT ─────────────────────────────────────────────────────
   Mirrors SUBROUTINE SGERG(X2, X3, HS, RM, X5, P, TC, Z, D)

   Inputs
     X3   CO2 mole fraction          (0.0 – 0.30)
     HS   calorific value  [MJ/m³]   (20  – 48,  metering 0 °C / 1.01325 bar)
     RM   relative density           (0.55 – 0.90, relative to air)
     X5   H2  mole fraction          (0.0 – 0.10)
     P    pressure  [bar]            (0   – 120)
     TC   temperature  [°C]          (-23 – 65)

   Returns { Z, D [mol/L = kmol/m³], X2 [calculated N2 mole fraction] }    */
function sgerg(X3, HS, RM, X5, P, TC) {
    if (P  <  0.0 || P  > 120.0) throw new Error('Pressure out of range (0 – 120 bar)');
    if (TC < -23.0 || TC >  65.0) throw new Error('Temperature out of range (−23 – 65 °C)');
    return SGERG1(P, TC, X3, X5, HS, RM);
}

window.sgerg = sgerg;
