# Math constants in JavaScript

List with 190+ math constants.

## Preview

```ts
/**
 * Antiquity
 */
const π = 3.14159265358979323846; // Pi, π
const PI = 3.14159265358979323846; // Pi, π
const τ = 6.283185307179586476925; // Tau, τ
const TAU = 6.283185307179586476925; // Tau, τ
const SQRT2 = 1.4142135623730950488; // Square root of 2, Pythagoras constant, √2
const SQRT3 = 1.73205080756887729352; // Square root of 3, Theodorus constant, √3
const GOLDENRATIO = 1.6180339887498948482; // Phi, The Golden Ratio, ϕ
const SQRT5 = 2.2360679774997896964; // Square root of 5, √5
const DELIAN = 1.25992104989487316476; // Cube root of 2, Delian Constant, 3√2

/**
 * Medieval and Early Modern
 */
const BUFFON = 0.63661977236758134307; // Buffon constant, 2/π
const HERMITE_KEPLER = 0.74048048969306104116; // Hermite constant Sphere packing 3D Kepler conjecture, μ_K
const W = 2.09455148154232659148; // Wallis Constant, W
const E = 2.71828182845904523536; // Euler's number, e
const LN2 = 0.69314718055994530941; // Natural logarithm of 2, ln(2)
const FAVARD_K1_WALLIS = 1.57079632679489661923; // Favard constant K1 Wallis product, π/2
const SQRT2PI = 2.50662827463100050241; // Square root of Pi times two, √2π
const SOPHOMORE1 = 0.78343051071213440705; // Sophomore's dream_1, I_1
const SOPHOMORE2 = 1.2912859970626635404; // Sophomore's dream_2, I_2
const LEMNISCATE = 2.62205755429211981046; // Lemniscate constant, varpi
const EULER_MASCHERONI = 0.5772156649015328606; // Euler–Mascheroni constant, γ (gamma)
const I_TO_I = 0.20787957635076190854; // i^ì, i to the power of i
const EB = 1.60669515241529176378; // Erdős–Borwein constant, E_B
const ET = 1.94359643682075920505; // Euler Totient constant, ET
const LAPLACE = 0.66274341934918158097; // Laplace limit, λ (lambda)
const GAUSS = 0.83462684167407318628; // Gauss's constant, G
const PISQUARED = PI * PI; // Pi squared, π^2
const FOURTHROOTFIVE = 1.49534878122122054191; // Fourth root of five, 4√5
const JOHN = 4.81047738096535165547; // John constant, γ (gamma), 1/(i^i)

/**
 * 19th century
 */
const R2 = 1.56155281280883027491; // The Triangular root of 2, R_2
const RAMANUJAN_SOLDNER = 1.45136923488338105028; // Ramanujan–Soldner constant, μ (mu)
const HERMITE = 1.15470053837925152901; // Hermite constant,  γ_2
const RIEMANN_ZETA = 1.64493406684822643647; // Riemann Function Zeta(2), ζ(2)
const LIOUVILLE = 0.110001000000000000000001; // Liouville number, £_Li
const R = 262537412640768743.999999999999250073; // Hermite–Ramanujan constant, R
const C = 0.91596559417721901505; // Catalan's constant, C
const DOTTIE = 0.73908513321516064165; // Dottie number, d
const M = 0.26149721284764278375; // Meissel-Mertens constant, M
const WEIERSTRASS = 0.47494937998792065033; // Weierstrass constant, σ (1/2)
const KASNER = 1.7579327566180045327; // Kasner number, R
const HAFNAR_SARNAK_MCCURLEY = 0.60792710185402662866; // Hafner–Sarnak–McCurley constant, 1/ζ(2)
const CAHEN = 0.64341054628833802618; // Cahen's constant, ξ_2
const P2 = 2.29558714939263807403; // Universal parabolic constant, P_2
const APERY = 1.20205690315959428539; // Apéry's constant, ζ(3)
const GELFOND = 23.1406926327792690057; //Gelfond's constant, e^π

/**
 * 1900–1949
 */
const L1 = 1.43599112417691743235; // Lebesgue constant (interpolation), L_1
const FAVARD = 1.23370055013616982735; // Favard constant, 3/4 * ζ(2)
const GOLDENANGLE = 2.39996322972865332223; // Golden angle, b
const K = 2.58498175957925321706; // Sierpiński's constant, K
const NIELSEN_RAMANUJAN = 0.82246703342411321823; // Nielsen–Ramanujan constant, ζ(2)/2
const L2 = 1.64218843522212113687; // Lebesgue constant L2, L2
const MANDELBROT_AREA = 1.5065918849; // Area of the Mandelbrot fractal, γ (gamma)
const GIESEKING = 1.01494160640965362502; // Gieseking constant,  π ln β
const BERNSTEIN = 0.28016949902386913303; // Bernstein's constant, β (beta)
const TR = 0.98770039073605346013; // Area bounded by the eccentric rotation of Reuleaux triangle, T_R
const TWIN_PRIMES = 0.66016181584686957392; // Twin Primes Constant, C_2
const SPIRAL_OF_THEODORUS = 1.86002507922119030718; // Spiral of Theodorus, constant of Theodorus, ∂ (partial)
const PLASTIC = 1.32471795724474602596; // Plastic number, ῥ (rho)
const L = 0.54325896534297670695; // Bloch-Landau constant, L
const GOLOMB_DICKMAN = 0.62432998854355087099; // Golomb–Dickman constant, λ (lambda)
const CFT = 0.66131704946962233528; // Feller–Tornier constant
const C10 = 0.12345678910111213141; // Base 10 Champernowne constant, C_10, C_FT
const GGS = 2.66514414269022518865; //Gelfond–Schneider constant, G_GS
const K0 = 2.6854520010653064453; // Khinchin's constant, K_0
const KHINCHIN_LEVI1 = 1.18656911041562545282; // Khinchin–Lévy constant, γ (gamma)
const KHINCHIN_LEVI2 = 3.27582291872181115978; // Khinchin–Lévy constant, γ (gamma)
const LEVY2 = 2.37313822083125090564; // Lévy 2 constant, 2 ln γ (gamma)
const C5 = 0.82699334313268807426; // Disk Covering, C_5
const CCR = 1.55138752454832039226; // Calabi triangle constant, C_CR
const MILLS = 1.30637788386308069046; // Mills' constant, θ (theta)
const EULER_GOMPERTZ = 0.59634736232319407434; // Euler–Gompertz constant, G

/**
 * 1950-1999
 */
const VAN_DER_PAUW = 4.53236014182719380962; // Van der Pauw constant, α (alpha)
const BETA = 1.38135644451849779337; // Beta, Kneser-Mahler polynomial constant, β (beta)
const LOCHS = 0.97027011439203392574; // Lochs constant £_Lo
const TANH1 = 0.76159415595576488811; //Hyperbolic tangent of 1, tanh 1
const W2D = 1.53960071783900203869; // Lieb's square ice constant, W_2D
const PI_TO_E = 22.45915771836104547342; // π^e, Pi to e
const NIVEN = 1.70521114010536776428; // Niven's constant, C
const BAKER = 0.83564884826472105333; // Baker constant, β_3
const BAXTER_FOUR_COLORING = 1.46099848620631835815; // Baxter's Four-coloring constant, C^2
const GAUSS_KUZMIN_WIRSING = 0.30366300289873265859; // Gauss–Kuzmin–Wirsing constant, λ_2
const PORTER = 1.46707807943397547289; // Porter's constant, C
const FEIGENBAUM_DELTA = 4.66920160910299067185; // Feigenbaum constant δ, δ (delta)
const CHAITIN = 0.0078749969978123844; // Chaitin's constants, Ω (omega)
const INFINITE_PRODUCT_AG = 1.75874362795118482469; // Infinite product constant, with Alladi-Grinstead, P_r1
const ALLADI_GRINSTEAD = 0.80939402054063913071; // Alladi–Grinstead constant, A_AG
const F = 2.80777024202851936522; // Fransén–Robinson constant, F
const ROBBINS = 0.66170718226717623515; // Robbins constant, Δ(3)
const FEIGENBAUM_ALPHA = 2.50290787509589282228; // Feigenbaum constant α, α (alpha)
const FRACTAL_DIMENSION_CANTOR = 0.63092975357145743709; // Fractal dimension of the Cantor set, d_f(k)
const CONNECTIVE = 1.84775906502257351225; // Connective constant, μ (mu)
const SALEM = 1.17628081825991750654; // Salem number, Lehmer's conjecture, σ_10
const PL = 1.96285817320964578286; // Reciprocal Lucas constant, P_L
const CHEBYSHEV = 0.59017029950804811302; // Chebyshev constant, λ_Ch
const CONWAY = 1.30357726903429639125; // Conway constant, λ (lambda)
const PREVOST = 3.35988566624317755317; // Prévost constant Reciprocal Fibonacci constant, ψ (psi)
const B2 = 1.902160583104; // Brun 2 constant = Σ inverse of Twin primes, B_2
const VARDI = 1.26408473530530111307; // Vardi constant, V_c
const Q = 0.28878809508660242127; // Flajolet and Richmond, Q
const HAFNER_SARNAK_MCCURLEY = 0.35323637185499598454; // Hafner–Sarnak–McCurley constant, σ (sigma)
const FRACTAL_DIMENSION_APOLLONIAN_TD = 1.305686729; // Fractal dimension of the Apollonian packing of circles by Thomas & Dhar, ɛ (varepsilon)
const FRACTAL_DIMENSION_APOLLONIAN_MULLEN = 1.305688; // Fractal dimension of the Apollonian packing of circles by McMullen, ɛ (varepsilon)
const B = 1.45607494858268967139; // Backhouse's constant, B
const GROTEHNDIECK = 1.78221397819136911177; // Gorthendieck constant, K_R
const MURATA = 2.82641999706759157554; // Murata Constant, C_m
const S1 = 1.09317045919549089396; // Smarandache Constant 1ª, S_1
const VISWANATH = 1.1319882487943; // Viswanath constant, C_Vi
const TIME = 0.6321205588285576784; // Time constant, tau, τ
const KOMORNIK_LORETI = 1.78723165018296593301; // Komornik-Loreti constant, q
const K_HARMONIC = 1.74540566240734686349; // Khinchin harmonic mean, K_-1
const REGULAR_PAPERFOLDING = 0.85073618820186726036; // Regular paperfolding sequence, P_f
const MADELUNG = 5.97798681217834912266; // Madelung Constant 2, C_2(2)
const ARTIN = 0.37395581361920228805; // Artin constant, C_Artin
const MRB = 0.18785964246206712024; // MRB constant, C_MRB
const DIMER_2D = 0.29156090403081878013; // Dimer constant 2D, Domino tiling, C/π
const HALL_MONTGOMERY = 0.17150049314153606586; // Hall-Montgomery Constant, δ_0
const SOMOS_QUADRATIC = 1.66168794963359412129; //Somos' quadratic recurrence constant, σ (sigma)
const STEINER = 1.44466786100976613365; // Steiner number, Iterated exponential Constant, e√e
const PLOUFFE_A = 0.15915494309189533576; // Plouffe's A constant, A, 1/2π
const GAUSS_LEMNISCATE = 1.85407467730137191843; // Gauss' Lemniscate constant, L/√2

/**
 * 2000-today
 */
const FOIAS_ALPHA = 1.18745235112650105459; // Foias constant_α (alpha), F_α
const FOIAS_BETA = 2.2931662874118610315; // Foias constant_β (beta), F_β
const GOH_SCHMUTZ = 1.11786415118994497314; // Goh-Schmutz constant, C_GS
const PELL = 0.58057755820489240229; // Pell constant, P_Pell
const CAREFREE = 0.70444220099916559273; // Carefree constant_2, C_2
const FIBONACCI_FACTORIAL = 1.22674201072035324441; // Fibonacci Factorial constant, F
const DOUBLE_FACTORIAL = 3.05940740534257614453; // Double factorial constant, C_n!!
const VOLUME_REULEAUX = 0.42215773311582662702; // Volume of Reuleaux tetrahedron, V_R
const RAABE = 0.91893853320467274178; // Raabe's formula, ζ'(0)
const JJGJJG = 0.9288358271; // Sum of the reciprocals of the averages of the twin prime pairs, JJGJJG, B_1
const SILVERMAN = 1.78657645936592246345; // Silverman constant, S_m
const KEPLER_BOUWKAMP = 0.1149420448532962007; // Kepler–Bouwkamp constant, ῥ (rho)
const PROUHET_THUE_MORSE = 0.41245403364010759778; // Prouhet–Thue–Morse constant, τ (tau)
const GAMMA_3_4 = 1.22541670246517764512; // Gamma(3/4), tfrac(3/4)
const HEATH_BROWN_MOROZ = 0.0013176411548531781; // Heath-Brown–Moroz constant, C_HBM
const KEMPNER_SERIE = 23.103447909420541616; // Kempner Serie(0), K_0
const C1 = 0.98943127383114695174; // Lebesgue constant, C_1
const C2 = 0.19452804946532511361; // 2nd du Bois-Reymond constant, C_2
const LUROTH = 0.78853056591150896106; // Lüroth constant, C_L
const STEPHENS = 0.57595996889294543964; // C_S
const TANIGUCHI = 0.67823449191739197803; // Taniguchi constant, C_T
const SILVER_ROOT = 3.24697960371746706105; // Silver root Tutte–Beraha constant, ς (varsigma)
const AREA_FORD = 0.87228404106562797617; // Area of Ford circle, A_CF
const CUBIC_RECURRENCE = 1.15636268433226971685; // Cubic recurrence constant, σ_3 (sigma 3)
const FRODA = 6.58088599101792097085; // Froda constant, 2^e
const MASSER_GRAMAIN = 0.64624543989481330426; // Masser-Gramain constant, C
const GIBBS = 1.85193705198246617036; // Gibbs constant, Si(π), Sin integral
const COPELAND_ERDOS = 0.23571113171923293137; // Copeland–Erdős constant, C_CE
const HAUSDORFF = 1.58496250072115618145; // Hausdorff dimension / Sierpinski triangle, log_2(3)
const V_8 = 2.02988321281930725004; // Figure eight knot hyperbolic volume, V_8
const CARLSON_LEVIN = 1.77245385090551602729; // Carlson–Levin constant
const STRONGLY_CAREFREE = 0.2867474284344787341; // Strongly Carefree constant, K_2
const MAGIC_ANGLE = 0.955316618124509278163; // Magic angle, θ_m (theta_m)
const RANDOM_WALK = 0.34053732955099914282; // Pólya Random walk constant, p(3)
const LANDAU_RAMANUJAN = 0.76422365358922066299; // Landau–Ramanujan constant, K
const UPPER_ITERATED = 0.69034712611496431946; // Upper iterated exponential, H_2n+1
const LOWER_ITERATED = 0.6583655992; // Lower límit iterated exponential, H_2n
const WRIGHT = 1.92878; // Wright constant, omega
const B4 = 0.870588379975; // Brun_4 constant = Σ inv.prime quadruplets, B_4
const R5 = 2.74723827493230433305; // Ramanujan nested radical, R5

/**
 * Other constants
 */
const BRONZE_RATIO = 3.30277563773199464655; // Bronze ratio, σ_Rr (sigma_Rr)
const SQRTTAUE = 4.13273135412249293846; // Square root of Tau*e, √τe
const TETRANACCI = 1.92756197548292530426; // Tetranacci constant, τ (tau)
const DAVICCI = 1.00743475688427937609; // DeVicci's tesseract constant, f_(3,4)
const TETRATION_I = 0.56755516330695782538; // Module of Infinite Tetration of i
const PRIMORDIAL = 0.70523017179180096514; // Primorial constant, Sum of the product of inverse of primes, P_#
const PLOUFFE_GAMMA = 0.14758361765043327417; // Plouffe's gamma constant, C
const SARNAK = 0.7236484022982000094; // Sarnak constant, C_sa
const BARNES_G = 1.78107241799019798523; // Exp.gamma, Barnes G-function, e^γ (e to gamma)
const RENYI_PARKING = 0.74759792025341143517; // Rényi's Parking Constant, m
const PARIS = 1.09864196439415648573; // Paris Constant, C_Pa
const EXPONENTIAL_FACTORIAL = 1.61111492580837673611; // Exponential factorial constant, S_Ef
const GLAISHER_KINKELIN = 1.28242712910062263687; // Glaisher–Kinkelin constant, A
const TROTT = 0.10841015122311136151; // Trott constant, T_1
const FRACTAL_DIMENSION_BOUNDARY_DRAGON_CURVE = 1.52362708620249210627; // Fractal dimension of the boundary of the dragon curve, C_d
const TRIBONACCI = 1.83928675521416113255; // Tribonacci constant, ϕ_3 (phi_3)
const GUMBEL_MEDIAN = 0.36651292058166432701; // Median of the Gumbel distribution, ll_2
const PI_TO_PI = 36.4621596072079117709; // Pi^pi, π^π
const GAMMA_1_4 = 3.62560990822190831193; // Gamma(1/4)
const GAMMA_1_6 = 5.56631600178023520425; // Gamma(1/6)
const ZETA_4 = 1.08232323371113819151; // Zeta 4, ζ(4)
const ZETA_5 = 1.03692775514336992633; // Zeta 5, ζ(5)
const IOACHIMESCU = 0.53964549119041318711; // Ioachimescu constant, 2 + ζ(1/2) (2 + zeta(1/2))
const FRACTAL_DIMENSION_KOCH_SNOWFLAKE = 1.26185950714291487419; // Fractal dimension of the Koch snowflake, C_k
const GOLDENSPIRAL = 1.3584562741829884352; // Golden spiral, c
const S5 = 2.79128784747792000329; // Nested radical S_5
const BESSEL = 0.697774657964007982; // Continued fraction constant, Bessel function, C_CF
const SUPER_SQRT2 = 1.55961046946236934997; // Super-Square root of 2, √2_s
const SUPER_SQRT3 = 1.82545502292483004004; // Super-Square root of 3, √3_s
const SUPER_CBRT2 = 1.47668433735786994708; // Super-Cube root of 2, 3√2_s
const OMEGA = 0.56714329040978387299; // Omega constant, Lambert W function, Ω (omega)
const FEIGENBAUM_DELTA3 = 5.96796870377745104099; // Feigenbaum constant δ of order 3, δ(3)
const E_TO_E = 5.154262241479259; // e to the power of e, e^e
```
