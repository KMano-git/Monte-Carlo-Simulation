# D+D+ elastic-collision fit memo for Codex

このメモは、D+D+ 弾性衝突の実装用に、論文中の**数式形**と**係数**をまとめたものです。
主に以下を含みます。

1. Krstić & Schultz (IAEA Volume 8) の D+ + D データ
   - 積分断面積フィット式
   - D+ + D 積分断面積フィット係数
   - D+ + D 微分断面積フィット式
   - D+ + D 微分断面積係数（raw extract）
2. Bachmann & Reiter (1995) の H+ + H 古典フィット
   - スケーリングの考え方
   - H+ + H cross-section fit coefficients (raw extract)

---

## 1. Krstić: 積分断面積の基本形

\[
\sigma_{\mathrm{el},\mathrm{mt},\mathrm{vi},\mathrm{se}}(E)
=
\frac{\sum_i a_i (\ln E)^i}{1+\sum_j b_j (\ln E)^j}
\]

- \(E\): center-of-mass collision energy [eV]
- 出力断面積単位: atomic unit (a.u.)
- 1 a.u. \(= a_0^2 = 2.80028\times 10^{-17}\,\mathrm{cm^2}\)

### D+ + D integral-fit coefficients (clean transcription)

#### Elastic
- a0-a3 = 0.637254E+03, -0.125310E+03, 0.167078E+02, -0.990209E+00

#### Momentum Transfer
- a0-a4 = 0.347935E+03, 0.419075E+03, 0.117426E+03, -0.324821E+01, -0.105819E+01
- b1-b3 = 0.131735E+01, 0.482389E+00, 0.401043E-01

#### Viscosity
- a0-a5 = 0.931543E+02, -0.391466E+02, 0.860210E+01, -0.595316E+01, 0.598129E+00
- b1-b4 = -0.139543E+00, 0.164132E+00, 0.334686E+00, 0.142442E+02
- b5 = 0.372665E-01

#### Spin Exchange
- a0-a5 = 0.173916E+03, -0.176678E+02, 0.507140E+00, 0.237852E+00, -0.234413E-01, -0.679823E+00

(注) 上の clean transcription は PDF の layout を見て再整形したもの。Momentum-transfer と viscosity は段組み配置のため raw extract も下に残す。

---

## 2. Krstić: 微分断面積の基本形

\[
2\pi \sin\theta \frac{d\sigma_{\mathrm{el}}}{d\Omega}(\theta)
=
\left[A+B(1-\cos\theta)+C\sin^2\theta\right]
\exp\left[
\frac{\sum_i a_i (\ln\theta)^i}{1+\sum_j b_j (\ln\theta)^j}
\right]
\]

- \(\theta\): CM scattering angle [rad]
- 出力単位: a.u. \(= a_0^2\,\mathrm{srad}^{-1}\)
- spin-exchange differential cross section では \(B=C=0\)

D+ + D の differential-fit coefficients はエネルギーごとに別係数です。
raw extract をそのまま下に添付します。Codex でパースして使う想定です。

---

## 3. Bachmann: スケーリングと cross-section fit

### 3.1 同位体スケーリングの考え方
Bachmann は、**中心質量系で interaction velocity が同じなら cross section は isotope-independent** とする単純質量スケーリングを採用している。

### 3.2 H+ + H cross-section fit
中間領域では

\[
\ln \sigma^{(l)}(E_r) \approx \sum_{n=0}^{8} a_n (\ln E_r)^n
\]

という多項式形を採用。
Appendix 2 に H+ + H の係数表あり。
raw extract を下に添付する。

---

# Raw extracted pages


## Krstic PDF page 104

1.2 
D+ + D —> D+ + D 
Energy (CM) 
(eV) 
0.1000 
0.1995 
0.5012 
1.0000 
1.9950 
5.0120 
10.0000 
19.9500 
50.1200 
100.0000 
Elastic 
(a.u.) 
9.773961E+02 
9.096440E+02 
7.024936E+02 
6.461367E+02 
5.621003E+02 
4.986651E+02 
4.488888E+02 
3.923312E+02 
3.492420E+02 
3.011882E+02 
Analytic fitting function 
Cross Section 
Momentum Transfer 
(a.u.) 
4.120358E+02 
4.008971E+02 
3.630347E+02 
3.487627E+02 
3.241112E+02 
2.897988E+02 
2.654433E+02 
2.419854E+02 
2.121037E+02 
1.907672E+02 
Viscosity 
(a.u.) 
1.813288E+02 
1.358733E+02 
1.079634E+02 
9.545577E+01 
6.752089E+01 
2.707555E+01 
1.410217E+01 
7.248538E+00 
2.436365E+00 
8.959529E-01 
Spin Exchange 
(a.u.) 
2.213106E+02 
2.056762E+02 
1.802378E+02 
1.752253E+02 
1.617602E+02 
1.448380E+02 
1.326478E+02 
1.209559E+02 
1.060618E+02 
9.538493E+01 
Vel,mt,m,se(E) 
= 
fc>(/n(£))A 
/ ( 1. + ] T b , ( / n ( £ ) ) M 
a.U, 
where E is the collision energy in the center of mass (CM) system expressed in eV and the 
cross section is in atomic units (1 a.u. = a\ = 2.80028E-17 cm2) 
Fitting parameters 
Elastic 
ao-a3: 
.637254E+03 
-.125310E+03 
.167078E+02 
-.990209E+00 
Momentum Transfer 
ao-a3: 
.347935E+03 
.419075E+03 
.117426E+03 
a4: 
-.105819E+01 
b rb 3: 
.131735E+01 
.482389E+00 
.401043E-01 
-.324821E+01 
Viscosity 
ao-a3: 
.931543E+02 
-.391466E+02 
.860210E+01 
a4-a5: 
-.595316E+01 
.598129E+00 
b rb 4: 
-.139543E+00 
.164132E+00 
.334686E+00 
.142442E+02 
.372665E-01 
Spin Exchange 
ao-a3: 
.173916E+03 
-.176678E+02 
.507140E+00 
a4-a5: 
.237852E+00 
-.234413E-01 
-.679823E+00 
92 



## Krstic PDF page 105

D+ + D -> D+ + D 
3 
03 
O 
CD 
CO 
CO 
o 
o 
103« 
< 
m2 
101 
m° 
»_d 
*^*H 
^ S 
\^\ 
**&£ 
JHQ* 
' ^ n 
^ 
9-C 
3 c / r s 
" " ^ r -
^ * >! 
H*, 
K 
\ ^ 
&-« 
% 
^ 
N * ^ 
k 
>> 
Elastic 
Momentum Transfer 
Viscosity 
— — Spin Exchange 
o Computed from fit to D 
Computed from fit to T< 
I II 
I 
I 
"<*H-N 
Vc 
—v 
\ 
cs 
~s 
i i r i 
^e-
^ 
. 
> 
\ 
\ 
EN 
=>-
s. 
^ 
H^ 
K 
\ k 
tet\ 
\ r - j f 
-i 
10 
10 
10 
10' 
Energy in Center of Mass System (eV) 
93 



## Krstic PDF page 106

1.2 
D+ + D —• D+ + D 
Elastic and Spin Exchange Differential Cross Sections 
Analytic fitting function 
27rsin{6)^j^{0) 
=[A + B{1 - cos(O)) + Csin\9)\ 
exp 
dn 
I>(Zn(*))' / l. + 5X/n(0))' 
\i=0 
/ 
\ 
j=l 
a.u., 
where A. B, C', a,-, and bj are coefficients depending on the center of mass collision energy (E, 
eV) and scattering angle (9, radians) and the cross section is in atomic units (1 a.u. = a2
0 srad-1 
= 2.80028E-17 cm2 srad-1). Note that for the spin exchange (se) differential cross section, B 
and C are zero. 
Fitting parameters 
E = .1000 eV 
Elastic 
a-o~a-5-
bi-b 3: 
A,B,C: 
Spin Exchange 
ao-a4: 
b!-b 4: 
A: 
E = .1259 eV 
Elastic 
a0- a.$: 
bi-b 4: 
A,B,C: 
Spin Exchange 
3-0" &5 • 
b!-b 4: 
A: 
E = .1585 eV 
Elastic 
ao-a5: 
bi-b 3: 
A,B,C: 
Spin Exchange 
ao-as'-
bj-b 4: 
A: 
.467625E+01 
-.320249E+00 
.102214E+01 
.370895E+01 
-.264780E+00 
.112631E+01. 
.441911E+01 
-.423590E+00 
.103422E+01 
.368426E+01 
-.539876E+00 
.117438E+01 
.465265E+01 
-.388821E+00 
.105376E+01 
.354822E+01 
-.368145E+00 
.115288E+01 
-.154536E+01 
-.180574E+00 
.550177E-01 
-.351569E+00 
-.240352E+00 
-.167191E+01 
-.223744E+00 
.202813E+00 
-.177531E+01 
-:288025E+00 
-.183310E+01 
-.214909E+00 
.190445E+00 
-.813437E+00 
-.240886E+00 
-.562238E+00 
-.449252E-01 
-.114432E+00 
-.929971E+00 
-.484357E-01 
-.794796E+00 
-.559973E-01 
-.312452E+00 
-.132486E+01 
.879424E-03 
-.791552E+00 
-.474349E-01 
-.460127E+00 
-.884318E+00 
-.398285E-01 
-.341237E+00 
-.328939E+00 
-.872631E-03 
-.485301E+00 
.234539E-02 
.377422E-01 
.139590E-01 
-.353663E+00 
-.240323E+00 
-.800510E-04 
-.441537E-02 
-.278597E-01 
.222773E-01 
.137280E+00 
-.208239E-02 
-. 119994 E-01 
.118090E-02 
.358255E-02 
.148016E-01 
.131805E-02 
.797220E-03 
94 



## Krstic PDF page 107

E = .1995 eV 
Elastic 
ao-a^: 
bi-b3: 
A, B, C: 
Spin Exchange 
ao-a-j: 
bi-b4: 
A: 
E = .2512 eV 
Elastic 
ao-as: 
bi-b4: 
A,B,C: 
Spin Exchange 
ao-as: 
bi-b4: 
•A: 
.443008E+01 
-.264761E-01 
.103866E+01 
.345519E+01 
-.401677E+00 
.121483E+01 
.417954E+01 
-.391009E+00 
.101381E+01 
.351407E+01 
-.493496E+00 
.106344E+01 
.168453E+00 
-.166518E-01 
.460375E-01 
•.909199E+00 
•.263802E+00 
.154441E+01 
.186822E+00 
.306258E-01 
.266725E+01 
.203150E+00 
.418292E+00 
-.226822E-01 
-A20366E+00 
-968169E+00 
-.431264E-01 
.393896E+00 
-.454280E-01 
-.616720E-01 
.243135E+00 
-.181665E-01 
.259378E+00 
-.969458E-03 
.196713E-01 
-.330737E+00 
.130829E-01 
.118582E-02 
.233075E-02 
.776751E+00 
.439025E+00 
.207473E+00 
.187202E-01 
.487794E-01 
.179713E-01 
E = .3162 eV 
Elastic 
ao-a2: 
bi-b5: 
A,B,C: 
Spin Exchange 
ao-ai': 
bi-b6: 
A: 
E = .3981 eV 
Elastic 
ao-as: 
bi-b4: 
A, B, C: 
Spin Exchange 
ao-ai: 
bi-b6: 
A: 
A22216E+01 
-.197936E+00 
.108099E+01 
.311341E+01 
-.658378E-01 
.102238E+01 
.415871E+01 
-.406091E+00 
.102941 E+Ol 
.298692E+01 
-.853280E-01 
.103413E+01 
.112520E+01 
-.102415E+00 
.164758E+00 
-.364947E-01 
-.545515E-01 
-.122528E-01 
.397117E+00 
-.158233E+00 
-.311893E-02 
-.926824E-04 
.464875E-01 
-.615046E-02 
-.396331E-03 
-.100617E-04 
-.194241E+01 
-.322310E+00 
-.392548E-01 
.500334E-01 
-.183028E+00 
-.274015E-01 
.226068E-02 
.153887E-01 
.118856E-02 
.350550E-02 
.365361E+00 
-.171858E+00 
-.474990E-01 
-.585560E-02 
-.351331E-03 
-.837274E-05 
E = .5012 eV 
Elastic 
ao-as: 
bi-b4: 
A,B,C: 
Spin Exchange 
bi-b4: 
A: 
E = .6310 eV 
Elastic 
ao-as: 
bi-b4: 
A,B,C: 
Spin Exchange 
ao-as: 
bi-b4: 
A: 
.420243E+01 
-.415115E+00 
.105164E+01 
.321243E+01 
-.382014E+00 
.121736E+01 
.421912E+01 
-.450216E+00 
.106584E+01 
.325277E+01 
-.572866E+00 
.108540E+01 
.233560E+01 
-.255364E+00 
.249909E+00 
.117928E+00 
.675454E-02 
.170110E+00 
-.775696E-02 
.553177E-02 
.281067E-01 
-.113961E+00 
.708100E+00 
-.854332E+00 
-.232093E+00 
-.175159E-01 
.257008E+00 
-.414328E-01 
-.985435E-03 
•.250655E+01 
-.454199E+00 
.311858E+00 
.163908E+00 
.943945E-02 
.183049E+00 
-.203906E-02 
.829227E-02 
.102316E+00 
-.280770E+00 
•.189788E+01 
-.658948E+00 
.446729E+00 
.180729E+00 
.155803E-01 
•.190543E+00 
.496967E-01 
.148566E-01 
95 



## Krstic PDF page 108

E = .7943 eV 
Elastic 
ao-as: 
bi-b4: 
A,B,C: 
Spin Exchange 
ao-as: 
bi-b4: 
A: 
E = 1.0000 eV 
Elastic 
ao-a5: 
bi-b4: 
A,B,C: 
Spin Exchange 
a0-a5: 
a6: 
bi-b5: 
A: 
.409783E+01 
-.426081E+00 
.106191E+01 
.295207E+01 
-.535692E+00 
.104693E+01 
.406031E+01 
-410627E+00 
.103802E+01 
.297629E+01 
-.730368E-01 
.149891E+01 
.103481E+01 
.255725E+01 
.166574E+00 
.706309E-01 
.175393E+01 
.159085E+00 
-.213183E+01 
-.187478E+00 
.644073E-01 
.455746E+01 
.349564E+00 
.222393E+00 
.467604E+00 
.174154E+00 
.946207E-02 
.613814E-02 
.822024E-02 
.194343E+00 
.261419E+00 
.516827E+00 
.173497E+00 
.136464E-01 
.438565E-01 
.128267E-01 
.377940E+00 
.204807E+00 
.109000E+00 
.623565E-02 
-.119303E-01 
.491171E-02 
•.157754E+00 
.166843E+01 
-.148422E+01 
-.191885E+01 
-.656324E+00 
-.107141E+01 
-.545901E+00 
-.729439E-01 
E = 1.2590 eV 
Elastic 
ao-as: 
bi-b4: 
A, B, C: 
Spin Exchange 
ao-a^: 
bi-b4: 
A: 
E = 1.5850 eV 
Elastic 
ao-as: 
bj-b 4: 
A,B,C: 
Spin Exchange 
ao-a2: 
bi-b4: 
A: 
.394037E+01 
-.415486E+00 
.101403E+01 
.173293E+01 
-.240274E+00 
.882905E+00 
.383428E+01 
-.477396E+00 
.950966E+00 
.138299E+01 
-.335400E+00 
.951405E+00 
-.197990E+01 
-.202787E+00 
.765877E-01 
.178299E+01 
-.542173E-01 
•.232640E+01 
•.202778E+00 
.875360E-01 
.159089E+01 
-.641218E-01 
.393652E+00 
.571923E-01 
.580143E-01 
.351054E-02 
-.240329E-01 
.203919E-02 
-.819521E-01 
.730684E-01 
-.406750E-02 
.110668E-03 
.353748E+00 
.133434E+00 
.864948E-01 
.503066E-02 
-.184402E-01 
.366969E-02 
-546438E-01 
-.528871E-02 
-.527104E-02 
-.156437E-03 
E = 1.9950 eV 
Elastic 
ao-as: 
t>i-b4: 
A, B, C: 
Spin Exchange 
ao-a2: 
bi-b4: 
A: 
E = 2.5120 eV 
Elastic 
ao-a5: 
b2-b4: 
A,B,C: 
Spin Exchange 
ao-ao: 
bi-b3: 
A: 
.366041E+01 
-.522065E+00 
.923605E+00 
.116326E+01 
-.612252E+00 
.108286E+01 
.342347E+01 
-.534692E+00 
.928147E+00 
.972234E+00 
-579903E+00 
.101079E+01 
-.275670E+01 
-.184909E+00 
.415083E-01 
.112642E+01 
-.840780E-02 
-.303967E+01 
-.164152E+00 
.205180E-01 
.122822E+01 
-.280724E-01 
-.150114E+00 
.390550E+00 
.150468E+00 
.811458E-02 
-.125630E-02 
.688250E-02 
.267389E+00 
-.511318E+00 
-.200295E-02 
-.147273E-03 
.141166E+00 
.651087E+00 
.197748E+00 
.101062E-01 
.141772E-01 
.888277E-02 
.362547E+00 
.356882E+00 
-.110977E-02 
96 



## Krstic PDF page 109

E = 3.1620 eV 
Elastic 
ao-as: 
a6: 
bi-b3: 
A,B,C: 
Spin Exchange 
ao-a2: 
bi-b3: 
A: 
E = 3.9810 eV 
Elastic 
&Q-SL5'-
bi-b4: 
A,B,C: 
Spin Exchange 
ao-a2: 
bi-b4: 
A: 
E = 5.0120 eV 
Elastic 
ao-as: 
bi-b4: 
A,B,C: 
Spin Exchange 
ao-as: 
bi-b4: 
A: 
E = 6.3100 eV 
Elastic 
ao-as: 
bi-b4: 
A,B,C: 
Spin Exchange 
ao-a2: 
bx-b4: 
A: 
E = 7.9430 eV 
Elastic 
ao-as: 
bi-b4: 
A,B,C: 
Spin Exchange 
ao-a2: 
bi-b2: 
A: 
E = 10.0000 eV 
Elastic 
ao-as: 
bi-b4: 
A,B,C: 
Spin Exchange 
ao-a2: 
bi-b4: 
A: 
.316889E+01 
.275336E-03 
-.379550E+00 
.942809E+00 
.891619E+00 
-.612784E+00 
.630038E+00 
.282169E+01 
-.491003E+00 
.952082E+00 
.557373E+00 
-.340863E+00 
.981279E+00 
.254215E+01 
-.471105E+00 
.986760E+00 
.558215E+00 
.583556E+00 
.884170E+00 
.228698E+01 
-.460872E+00 
.104659E+01 
.467819E+00 
-.558438E+00 
.101389E+01 
.228921E+01 
-.465475E+00 
.996528E+00 
.281764E-02 
-.740087E+00 
.131012E+01 
.203523E+01 
-.443507E+00 
.107318E+01 
-.436119E+00 
-.471215E+00 
.135804E+01 
-.243300E+01 
-'.277631E+00 
-.648701E-01 
.150894E+01 
-.524947E-02 
-.323451E+01 
-.987852E-01 
.478702E-01 
.192722E+01 
-.103266E+00 
-.306532E+01 
-.101752E+00 
.626973E-01 
.181315E+01 
-.118321E+01 
-.292619E+01 
-.107064E+00 
.808577E-01 
.112073E+01 
-.129359E+00 
-.293727E+01 
-.109195E+00 
.970789E-01 
.155187E+01 
-.600856E-02 
-.269731-E+01 
-.121372E+00 
.134148E+00 
.183175E+01 
-.127622E+00 
-.201571E+00 
-.376502E-01 
.494648E+00 
-.441995E+00 
.465892E-03 
.943567E+00 
.465729E-01 
.175790E+00 
.100477E+00 
-.109099E-01 
.112693E+01 
.393008E-01 
-.216795E-01 
.895246E+00 
-.450895E+00 
.125718E+01 
.319022E-01 
-.265551E+00 
-.489230E+00 
-.250963E-01 
.124523E+01 
.305704E-01 
-.591616E+00 
-.602923E+00 
.130512E+01 
.211236E-01 
-.800000E+00 
.104328E+00 
-.116121E-01 
.498195E+00 
.120195E+01 
.124707E-01 
-.407615E-03 
.115747E+01 
.988869E-02 
-.135841E+01 
.429270E+00 
v 
.109771E+01 
.759918E-02 
-.129491E-02 
.108302E+01 
.737070E-02 
.996913E+00 
.518345E-02 
-.373928E-03 
.143732E+00 
.287439E+00 
.247237E+00 
-.758284E+00 
.209849E+00 
.206125E+00 
.169495E+00 
.113705E-01 
.136415E-01 
.112969E-01 
.429247E+0C 
.920064E-02 
.896996E-02 
.703948E-02 
97 



## Krstic PDF page 110

E = 12.5900 eV 
Elastic 
ao-as: 
bi-b4: 
A,B,C: 
Spin Exchange 
ao-a2: 
bi-b6: 
A: 
E = 15.8500 eV 
Elastic 
ao-as: 
b r b 4 : 
A,B,C: 
Spin Exchange 
ao-ao: 
bi-b6: 
b7-bg: 
A: 
.203825E+01 
-.269264E+01 
.129558E+01 
.984729E+00 
.166761E+00 
-444716E+00 
-.122556E+00 
.202010E-01 
.503944E-02 
.985459E+00 
.169907E+00 
-.109834E+01 
.689291E-02 
-.574303E+00 
.142550E+01 
.262627E+00 
-309957E+00 
-.271542E+00 
-.619246E-01 
.102812E+01 
-.698839E-02 
-.372146E-03 
-.730809E-05 
Warning: Fitted elastic differential cross section does not accurately yield o~vi 
.204447E+01 
-.268261E+01 
.126759E+01 
.961271E+00 
.161377E+00 
.661141E-02 
-.446287E+00 
-.126933E+00 
.181179E-01 
.473316E-02 
.942656E+00 
.178649E+00 
-.100000E+01 
-.683210E+00 
.191221E+01 
.424795E+00 
.217575E+00 
-.209987E+00 
-.276587E+00 
-.129063E+00 
-.304013E-01 
-.405727E-02 
-.311921E-03 
-.129069E-04 
-.223016E-06 
.117802E+01 
E = 19.9500 eV 
Elastic 
ao-as: 
bi-b4: 
A,B,C: 
Spin Exchange 
ao-a4: 
bi-b6: 
A: 
E = 25.1200 eV 
Elastic 
ao-as: 
bi-b4: 
A,B,C: 
Spin Exchange 
ao-a4: 
bi-b5: 
A: 
Warning: Fitted elastic differentia! cross section does not accurately yield crv{ 
.204513E+01 
-.268325E+01 
.126992E+01 
-.447414E+00 
-.126187E+00 
.179698E-01 
.861309E+00 
.208200E+00 
-.100000E+01 
.957310E+00 
.160519E+00 
.469870E-02 
-.109114E+01 
.103367E+01 
.124220E+01 
.587868E+00 
.254683E+01 
.111243E+01 
-.809625E+00 
-.455694E+00 
-.777382E+00 
.521296E-02 
.354899E-03 
Warning: Fitted elastic differential cross section does not accurately yield avi 
155529E+00 
.205085E+01 
-.267280E+01 
.124566E+01 
-.448688E+00 
-.129883E+00 
.161366E-01 
.820624E+00 
.213632E+00 
-.100000E+01 
-.144405E+01 
.123458E+01 
.921973E-01 
.368104E+00 
.117371E+00 
-.780350E+00 
.128220E+01 
.935796E+00 
.442277E-02 
.287893E+01 
-.814198E+00 
.246119E-02 
.675589E-04 
.655543E-02 
.898222E-05 
.629641E-02 
E = 31.6200 eV 
Elastic 
ao-as: 
bi-b4: 
A,B,C: 
Spin Exchange 
ao-a^: 
bi-b4: 
A: 
E = 39.8100 eV 
Elastic 
ao-as-. 
b]-b6: 
A,B,C: 
Spin Exchange 
ao-a2: 
bi-b3: 
A: 
Warning: Fitted elastic differential cross section does not accurately yield o\i 
.204942E+01 
-.268024E+01 
.125797E+01 
.937258E+00 
.155826E+00 
-.450513E+00 
-.127576E+00 
.165056E-01 
.445455E-02 
.756992E+00 
.231059E+00 
-.900000E+00 
.628381E-02 
-.224695E+01 
.260816E+01 
.414866E+00 
-.200755E-01 
-.500666E+00 
-.181194E+00 
-.163103E-01 
-.367304E-03 
.112343E+01 
.397848E-02 
.247437E+00 
-.175729E+01 
.115825E+01 
.761900E+00 
.115393E+00 
.441212E-02 
•.319955E+00 
-.290263E+00 
-.679068E-01 
-.985600E-02 
-.897343E-03 
-.225008E-04 
.904611E+00 
.495906E-01 
.263866E+00 
-.315069E+01 
.408915E+01 
-.103727E+00 
-.636083E+00 
-.509533E-01 
-.178664E-02 
.147048E+01 
98 



## Krstic PDF page 111

E = 50.1200 eV 
Elastic 
ao-as: 
bi-b6: 
A, B, C: 
Spin Exchange 
ag-a2: 
bi-b5: 
A: 
E- 
63.1000 eV 
Elastic 
ao-as: 
bi-b6: 
A,B,C: 
Spin Exchange 
ao-a2: 
bi-b4: 
A: 
.248984E+00 
-.327135E+00 
.850868E+00 
-.306105E+01 
-.616720E+00 
.876634E+00 
.250968E+00 
-.349475E+00 
.775014E+00 
-.411091E+01 
-.633308E+00 
.743076E+00 
-.176247E+01 
-.287696E+00 
.643532E-01 
.315806E+01 
-.129115E+00 
-.179219E+01 
-.279897E+00 
.100100E+00 
.517558E+01 
-.809774E-01 
.116661E+01 
.756407E+00 
.113242E+00 
.426666E-02 
-.657595E-01 
-.933116E-02 
-.836493E-03 
-.202568E-04 
.320403E+00 
.340093E+00 
-.270568E-02 
.664440E-03 
.355866E-04 
.119826E+01 
.756130E+00 
.110011E+00 
.396040E-02 
-.573517E-01 
-.719284E-02 
-.591235E-03 
-.961807E-05 
•.756506E+00 
•.469997E+00 
-.131608E-01 
-.597900E-03 
E = 79.4300 eV 
Elastic 
ao-as: 
bi-b5: 
A,B,C: 
Spin Exchange 
ao-a2.' 
b1-b4: 
A: 
E = 100.0000 eV 
Elastic 
ao-as: 
bi-b5: 
A,B,C: 
Spin Exchange 
ao-a2: 
bi-b4: 
A: 
Warning: Fitted elastic differential cross section does not accurately yield <jvi 
.259839E+00 
-.180394E+01 
.120472E+01 
.742606E+00 
.104606E+00 
.357125E-02 
-.367530E+00 
-.276366E+00 
-.514901E-01 
-.552782E-02 
-.379661E-03 
.708988E+00 
.105348E+00 
-.750000E+00 
-.434102E+01 
.397219E+01 
.501763E+00 
-.616269E+00 
-.138469E+00 
-.755848E-02 
-.135204E-03 
.922278E+00 
Warning: Fitted elastic differential cross section does not accurately yield avi 
.272160E+00 
-.178173E+01 
.117749E+01 
.715882E+00 
.990868E-01 
.330418E-02 
-.368500E+00 
-.279184E+00 
-.532186E-01 
-.579052E-02 
-.382924E-03 
.702837E+00 
.698164E-01 
-.750000E+00 
-.536273E+01 
.515131E+01 
.771888E+00 
-.491597E+00 
-.147314E+00 
-.118633E-01 
-.433425E-03 
.180480E+01 
99 



## Krstic PDF page 112

D+ + D -» D+ + D 
D+ + D -> D+ + D 
ECM = 0.1eV 
ECM = 0.1995 eV 
D 
"a 
<x> ioH 
• 
K 103 
CM 
102 
101 
10"5 10"4 10"3 10"2 10"1 10° 
1 0"5 10"4 10"3 10"2 10"1 
10° 
10u 
10 
-1 
10"' 
Fit to Elastic 
Elastic 
Fit to Spin Exchange 
Spin Exchange 
1 
' 
• 
0 
3 
0 
Scattering Angle in Center of Mass System (rad) 
100 



## Krstic PDF page 113

D+ + D -> D+ + D 
D+ + D -> D+ + D 
ECM = 0.5012 eV 
ECM = 1 eV 
a 
c 
w 1cr 
10" 1.0"* 10"° 10" 10" 10". 10~5 10"4 10*3 10"2 10"1 10° 
10< 
4 
i i i i i i i i I i I i i i i i I I I i i i I i i i i i I 
CM 
101 
10' 
.-1 
1 0 ' 
• 
10 
10"; 
-2 
10" 
10 
-5 
I I I I I I I I I II I I I I I I I I I I I I I I I I I I I II 
— Fit to Elastic 
— Elastic 
— Fit to Spin Exchange 
Spin Exchange 
* * * 
i . . • 
i .i 
0 
1 
2 
3 
0 
1 
2 
Scattering Angle in Center of Mass System (rad) 
101 



## Krstic PDF page 114

D+ + D -> D+ + D 
D+ + D -» D+ + D 
ECM = 1.995 eV 
ECM = 5.012 eV 
a 
•D 
<x> 
"co 
CM 
1 
10 10 10 10 10 10 10 10 10 10 10 10 
.1 I I U I I I I I I I I I I I I I I I I I I I I I I I I I i 
10 
10 
10 
-4 
• 
-6 
Fit to Elastic 
Elastic 
Fit to Spin Exchange 
Spin Exchange 
' • 
' 
* •* 
0 
1 
Scattering Angle in Center of Mass System (rad) 
102 



## Krstic PDF page 115

D+ + D -> D+ + D 
D+ + D -» D+ + D 
ECM=10eV 
a 10 
•6 
' 
i 11 I I ^ 
i niuq 
i 
i i iiiirf 
i mini 
" * • * 
" J 
- 1 
ECM = 19.95 eV 
• >*V 
. nft( 
I I 
r 
, 
" 
. 
* • 
S 
*s 
r 
s 
\ 
*f 
• 
^ 
4$S 
r 
A 
</ 
: 
^ 
</ 
" 
^ 
4S 
• 
, 
fix. 
* / r • 
// 
' 
^m*j 
• 
«^ 
JT 
" 
' 
-
* 
" 
1-1 
-5 
10 10 10 10 10 10 
10 10 10"° 10 10 10 
10 ki 
10" 
10' 
10" 
Fit to Elastic 
—— Elastic 
Fit to Spin Exchange 
Spin Exchange 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 
0 
1 
Scattering Angle in Center of Mass System (rad) 
103 



## Krstic PDF page 116

D+ + D -» D+ + D 
D+ + D -» D+ + D 
ECM = 50.12 eV 
EOM = 100eV 
§10-° 
-
1 
* • * 
" * 
-
1 
• • • • —
J 
• ' 
10"
5 10"4 10"3 10"2 10'1 10° 
10"5 10"4 10"3 10"2 10"1 10' 
1 0 " 
10"4 
ID'5 
10"6 
10"7 
10"8 
Fit to Elastic 
Elastic 
] 
Fit to Spin Exchange 
Spin Exchange 
• 1 1 1 1 1 1 1 1 1 1 1 
• 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ] 
J 
Spin Exchange 
| 
L 
J 
0 
1 
2 
3 
0 
1 
2 
3 
Scattering Angle in Center of Mass System (rad) 
104 



## Krstic PDF page 117

1. Hydrogen-ion-hydrogen-atom elastic collisions 
1.3 T + + T 
Important Note 
The calculations of both the differential and integral (elastic, momentum transfer, and 
viscosity) cross sections for the symmetric systems H + + H, D + + D, and T + + T have 
been performed assuming indistinguishability of the constituent nuclei. As a consequence, the 
elastic cross sections reported here are the coherent sum of the interfering processes of elastic 
scattering and spin exchange. Thus, this procedure results in a double counting if in a particular 
application of the data the spin exchange (resonant charge transfer) cross section needs to be 
treated separately. In that case the cross sections should be re-computed by subtracting the 
spin exchange differential cross section from the "total" elastic differential cross section and 
integrating this "pure" elastic differential cross section with the appropriate transport weighting 
functions to obtain the required moments. For example, the momentum transfer cross section 
tabulated here otherwise accounts for both spin exchange and elastic processes in the high 
energy limit. The error introduced by such a "decoupling" of the elastic and spin exchange 
cross sections does not exceed 10% in the energy range considered. (See the Introduction in 
Part A for more explanation.) 
105 


## Bachmann PDF page 52

96 
Contrib. Plasma Phys. 35 (1995) 1 
9 Appendix 2: Fit Coefficients for Cross sections 
Table 1 
Fitting parameterss for cross sections for R1 H+ + H 
SIGMA-TOTAL(ELAB) 
ELABMIN = ELABMAX = 1.8206eV 
a0 
0.000000000000D~00 
a1 
0.000000000000D+00 
a3 
0.000000000000D+00 
a4 
0.000000000000D+00 
a6 
0.000000000000D+00 
a7 
0.000000000000D+00 
a10 - .3253031352541D +01 
all 
-2.559032645641D -01 
arO - 3.2629373574MD+Ol 
arl - 8.719626183599D -02 
ELABMIN = 0.02 eV, ELABMAX = 200 eV 
a0 -3.349115100108D+01 
a1 -4.047040620920D- 11 
a3 -5.224890973622D-03 
a4 - 1.019115858754D-04 
a6 -4.33625901 1986D-05 
a7 - 1.781020734395D-06 
a10 
-3.320677627738D +01 
all 
-2.2059420401 12D -01 
arO -2.753878563969D +01 
arl -2.000000000000D +OO 
ELABMIN = 0.02 eV, ELABMAX = 200 eV 
SIGMA-DIFFUSION 
SIGM A-VISCOSITY 
a0 - 3.353420922048D +01 
a1 -3.522409780724D-01 
a3 -4.282561W6823D -03 
a4 - 3.230618998917D - 04 
a6 - 1.753965583282D-05 
a7 -4.580920664987D-07 
a10 -3.330015157525D+01 
all - 1.9926253664881)-01 
arO - 2.7093294272601) + 01 
ar 1 - 2.000000000000D + 00 
a2 
0.000000000000D+00 
a5 
0.000000000000D+00 
a8 
0.000000000000D+00 
a12 - 149996483552D -02 
ar2 -7.346647926269D-02 
a2 -4.340959073105D-02 
a5 - 3.3 141 57761518D -06 
a8 
1.2203935505627D-06 
a12 
0.000000000000D + 00 
ar2 
0.000000000000D + 00 
a2 - 3.587214262651D -02 
a5 -4.343173698940D-05 
a8 
3.738689325195D -07 
a12 
0.000000000000Df00 
ar2 
0.000000000000D+00 
Table 2 
Fitting parameters for cross sections for R2 H+ + He 
SIGMA-TOTAL(ELAB) 
ELABMIN = ELABC = 0.5081 eV, ELABMAX = ELABO = 29.4431 eV 
a0 -3.357907136508DfOl 
a1 -9.81 1659406594D-02 
a3 - 1.259671949006D+00 
a4 -4.4739475199848-02 
a6 - 1.203733922915D+00 
a7 
3.525830383820D-01 
a10 
-3.355838377904D+Ol 
all 
-2.845473342853D-01 
arO 
-3.706830076698D+Ol 
arl 
4.204258692619D-01 
ELABMIN = 0.0125 eV, ELABMAX = 125 eV 
a0 -3.4255853289531)+01 
a1 -8.999762959781D-01 
a3 
1.549750110754D-02 
a4 
3.963555202866D-02 
a6 -2.207534449376D-03 
a7 -3.378852519380D-05 
a10 - 3.390101844960D +01 
all 
-2.1 11706771 112D -01 
arO 
-3.034765152080D+Ol 
arl 
-2.000000000000D+00 
ELABMIN = 0.0125 eV, ELABMAX = 125 eV 
a0 - 3.443725345071D+Ol 
a1 -4.337427858507D-01 
a3 -6.451669335555D-02 
a4 
2.950009865269D-02 
a6 - 1.589840628629D - 03 
a7 - 1.502468439244D -04 
a10 
2.1 1 170677 1 1 12D - 01 
arO 
-2.978907423990D+Ol 
arl -2.000000000000D+W 
SIGMA-DIFFUSION 
SIGMA-VISCOSITY 
- 3.43227603 1579D + 01 
a1 1 
a2 
3.798308269292D-001 
a5 
1.565182597363D+00 
a8 - 3.668922671043D-02 
a12 - 1.351427675077D -02 
ar2 - 9.648359210100D - 02 
a2 - 3.434858 1248 1 1D - 01 
a5 
3.343570605088D-04 
a8 
4.224511209820D-05 
a12 
0.000000000000Df00 
ar2 
0.000000000000D + 00 
a2 -2.896488696126D-01 
a5 
5.752283385868D-03 
a8 
3.151 161681447D-05 
a12 
0.000000000000D+OO 
ar2 
0.000000000000D + 00 



## Bachmann PDF page 53

P. BACHMANN, 
D. REITER, Elastic Processes in Hydrogen-Helium Plasmas 
97 
Table 3 
Fitting parameters for cross sections for R3 H+ + H2 
SIGMA-TOTAL(ELAB) 
ELABMIN = ELABC = 1.5598 eV, ELABMAX = ELABO = 61.8164 eV 
a0 - 3.452141819446D +01 
a1 
1.092015526305D +01 
a3 
3.466297654768D + 01 
a4 - 2.524607958646D +01 
a6 -2.770065796605D +00 
a7 
3.796353200921D-01 
a10 - 3.275286840950D+01 
all 
-2.351764912137D-01 
arO - 3.537275807 146D + 01 
arl 
2.144573517210D - 01 
SIGMA-DIFFUSION 
ELABMIN = 0.015 eV, ELABMAX = 150 eV 
a0 -3.318680874597D+Ol 
a1 -3.580417289312D-01 
a3 -5.005702120342D-02 
a4 
2.369248748869D-02 
a6 - 1.357018742589D-93 
a7 - 1.39377609085533-04 
a10 
3.839304000000D- 15 
all 
-1.726918000000D-01 
arO 
2.5685026oooOOD - 12 
arl 
0.000000000000D + 00 
SIGM A-VISCOSITY 
ELABMIN = 0.015 eV, ELABMIN = 150 eV 
a0 -3.362402037774D+01 
a1 -2.337285826242D-01 
a3 -4.473235212373D-02 
a4 -4.691524784882D-03 
a6 
4.229065229431D-04 
a7 -6.739555319843D-05 
a10 
3.053906071449D- 15 
all - 1.726917299089D-01 
arO 
2.833830531873D- 12 
arl 
0.000000000000D+00 
a2 -2.732690257819D +01 
a5 
1.092376446349D + 01 
a8 -2.168988142310D-02 
a12 - 1.045602118569D-02 
ar2 -4.643079956637D-02 
a2 -2.27438237695 1D-01 
a5 
5.013459267775D -03 
a8 
3.029808591929D-05 
a12 
0.000000000000D+00 
ar2 
0.000000000000Dc00 
a2 - 5.404526201247D-02 
a5 
3.121568334037D-03 
a8 -7.756198335533D-06 
a12 
0.000000000000D+00 
ar2 
0.000000000000D+00 
Table 4 
Fitting parameters for cross sections for R4 He' + He 
SIGM A-TOTAL(ELAB) 
ELABMIN = ELABC = 1.2122eV, ELABMAX = ELABO = 64.6090eV 
a0 -3.336949020454D +01 
a1 
4.374909804779D +00 
a3 
2.345459194687D + 01 
a4 - 1.969436659467D + 01 
a6 -2.604153028956D+00 
a7 
3.801132783280D-01 
a10 - 3.2910713302481) +01 
all 
-2.416669402887D -01 
arO 
-3.6646919254241)+01 
arl 
4.752719886448D-01 
SIGMA-DIFUSION 
ELABMIN = 0.02 eV, ELABMAX = 200,eV 
a0 -3.332091557452D +01 
a1 -3.823354679977D-01 
a3 - 8.177418933677D -02 
a4 
2.593188019755D -02 
a6 - 1.649825718076D-03 
a7 -2.491587647454D-04 
a10 
4.524094200000D- 15 
all - 1.115060000000D-01 
arO 
7.6730207ooOOOD - 13 
arl 
0.000000000000D + 00 
SIGMA-VISCOSITY 
ELABMIN = 0.02 eV, ELABMAX = 200 eV 
a0 -3.379346231200D+Ol 
a1 - 1.740525006979D-01 
a3 -8,2238473151341)-02 
a4 -1.443276051210D-03 
a6 -5.5932944418448-005 
a7 - 1.742244159818D-04 
a10 
3.261848864837D- 15 
all - 1.115060177785D-01 
arO 
1.120563328874D-12 
arl 
0.000000000000D+00 
a2 - 1.517973301721D+01 
a5 
9.472303986781D + 00 
a8 -2.282922057203D-02 
a12 
-93213779217570 -03 
ar2 -8.280792916138D-02 
a2 -2.666453887008D-01 
a5 
8.320863897668D - 03 
a8 
4.351897658362D-05 
a12 
0.000000000000D+00 
ar2 
0.000000000000D + 00 
a2 -8.091712353563D-02 
a5 
6.530393601967D-03 
a8 
1.068285383642D - 05 
a12 
0.000000000000D+OO 
ar2 
0.000000000000D+00 



## Bachmann PDF page 54

98 
Table 5 
Contrib. Plasma Phys. 35 (1995) 1 
R(T) 
R ( T )  (cm**3/s) R1 H + H+ 
a0 - 1.820301019850D+ 1 
a3 -9.219441585651D-03 
a6 - 1.021209574651D-04 
R(T) (cm**3/s) R2 H + He 
a0 - 1.937488793273D+01 
a3 - 2.164290342772D - 02 
a6 
1.082584067859D-03 
R(T) (cm**3/s) R3 H2 + H+ 
a0 -1.852351664804D+Ol 
a3 - 1.354211329590D-02 
a6 
1.027842439916D -03 
R ( T )  (cm**3/s) R4 He + H, 
a0 - 1.927232246647D-01 
a3 - 5.6W616121964D - 03 
a6 
1.523437071234D -03 
a1 
2.129097882981D -01 
a4 
5.351258558196D-04 
a7 
1.044107275669D-05 
a1 
1.759 108971 539D - 01 
a4 - 1.3943958753511)-02 
a7 -2.740989886580D-04 
a1 
2.6977 18272304D - 0 1 
a4 -1.236691940121D-02 
a7 
-2.222408616640D-04 
a1 
1.490667894727D -01 
a4 - 1.908223060123D-02 
a7 - 3.1569909241 24D - 04 
a2 -3.319721546062D-02 
a5 
3.784499807334D -04 
a8 - 3.937241365267D - 07 
a2 - 1.693965011113D-02 
a5 
2.443778636181D-03 
a8' 
1.663374157798D-05 
a2 
1.350651671297D -02 
a5 
1.289363581693D -03 
a8 
1.252099405338D -05 
a2 
2.149231447853D -02 
a5 
1.394534352967D-03 
a8 
1.763329693714D-05 
Table 6 
OMEGA1 1 (T) 
OMEGAll (cm**3/s) R1 H + H+ 
a0 -2.067495970296D +01 
a1 -6.7544543866518-02 
a3 - 7.4299227822291) -03 
a4 - 7.808813003874D -05 
a6 -4.050666345184D -05 
a7 
1.836578913379D-05 
OMEGAll (cm**3/s) R2 H + He 
a0 -2.203507049378D +01 
a1 -6.666508829508D-01 
a3 
7.285987500237D -02 
a4 - 7.51 11370582OOD -03 
a6 
1.170265885717D -03 
a7 -7.290932193552D-05 
OMEGAll (cm**3/s) R3 H2 + H+ 
a0 -2.058948378358D+Ol 
a1 -3.091608300462D-01 
a3 
5.940618635367D -02 
a4 
9.580654599527D -03 
a6 - 2.415856443530D -04 
a7 
1.98193 1289049D - 04 
OMEGAll (cm**3/s) R4 He + H+ 
a0 -2.167965486890D+Ol 
a1 -5.460360301 179D-01 
a3 
6.762120367287D-02 
a4 
2.761214449593D -03 
a6 
3.881535032205D-04 
a7 
7.237252600152D-05 
a2 - 6.9061 861797 18D - 02 
a5 - 1.032158595103D-04 
a8 - 1.377982388469D -06 
a2 - 7.30163433276 1D - 02 
a5 -5.064240875152D-03 
a8 
3.61 1435091385D -07 
a2 - 1.157943457107D-01 
a5 - 5.46751 1099899D -03 
a8 - 1.443879934691D-05 
a2 - 1.181746519736D-01 
a5 - 5.300088703885D -03 
a8 -7.356271943058D-06 



## Bachmann PDF page 55

P. BACHMANN, 
D. BITER, 
Elastic Processes in Hydrogen-Helium Plasmas 
99 
Table 7 
OMEGA22(T) 
_ _ _ _ ~ ~ ~  
~~ 
~ 
_ _ _ _ ~  
OMEGA22 (cm++3/s) R1 H + H’ 
a0 - 1.970770010292D$Ol 
a1 -4.750547016654D-02 
a3 -5.020232904192D-03 
a4 - 1.897148429605D-03 
a6 
1.455 101056679D - 04 
a7 
6.1 137261 87722D -06 
OMEGA22 (cm++3/s) R2 H + He 
a0 - 2.11971 81 18212D + 01 
a1 - 6.047938325697D -01 
a3 
8.754381039703D-02 
a4 -4.745795404802D-03 
a6 
1.28694907373OD-03 
a7 -3.373371121283D-05 
OMEGA22 (cm*+3/s) R3 H, + H+ 
a0 - 1.9825201998890 +01 
a1 -2.045986618368D-01 
a3 
5.709736784682D -02 
a4 
1.535016658865D -02 
a6 -4.624966207800D -04 
a7 
2.9968933309641) -04 
OMEGA22 (cm++3/s) R4 He + He’ 
a0 - 2.086953660006D +01 
a1 -4.470227740985D -01 
a3 
7.173876575209D-02 
a4 
1.094500514641D-02 
a6 
2.665521010922D -05 
a7 
2.035371868694D-04 
a2 - 6.7 17043480643D - 02 
a5 - 5.978277546350D -04 
a8 - 1.495444636173D-06 
a2 - 1.105016966801D-01 
a5 - 7.12053844143 1 D -03 
a8 -3.452213060435D-06 
a2 - 1.398594524973D-01 
a5 - 7.1 133 16092217D -03 
a8 -2.193592594926D-05 
a2 - 1.654791308926D-01 
a5 -7.279575206851D-03 
a8 - 1.639359039345D-05 
References 
[l] REITER, D., in “Atomic and Plasma Material Interaction Processes”, eds. JANEV, R. K., DRAWN, 
121 SCHMIDER, 
R., BITER, 
D., ZEHRFELD, 
H. P., et al., J. Nucl. Mat. 1%-197 (1992) 810. 
[3] SCHNEIDER, 
R., BRAAMS, 
B., BITER, D., et al., Contrib. Plasma Phys. 32 (1992) 450. 
[4] KASTELEWICZ, 
H., SCHNEIDER, 
R., REITER, D., et al., Contrib. Plasma Phys. 32 (1992) 456. 
[5] WATKINS, 
M. L., REBUT, P. H., Proc. Inter. Conf. on Plasma Physics, 19th EPS, Innsbruck 1992, 
[6] BACHMANN, 
P., REITER, 
D., Max-Planck-Institut Report IPP 8/1 (May 1992). 
[8] REITER, D., BACHMANN, 
P., PRINJA, A. K., Contrib. Plasma Phys. 32 (1992) 261. 
[9] CUPINI, E., DE MA‘ITEIS, A., ENEA Report RT/TIB/87/33. 
H. W., Elsevier Sciences Publishers, 1993. 
16C, 11-731. 
[7] BACHMANN, P., REITER, 
D., PRINIA, A. K., J.  NU^. Mat. 196-198 (1992) 865. 
[lo] CUPINI, E., DE MATTEIS, 
A., IL NUOVO CIMENTO 11D (1089) 1489. 
[ll] ABOU-GABAL, 
H. H., EMMERT, 
G. A., Nucl. Fusion 31 (1991) 407. 
[12] HAAS, G., DUCHS, D., EHRENBERG, 
J., et al., EPS Berlin 1991, Contr. Papers 111-101; 
HAAS, G., BACHMANN, 
P., DUCHS, D., et al., J. Nucl. Mat. 196-198 (1992) 481. 
[I31 REITER, D., J. Nucl. Mat. 196-198 (1992) 80. 
[I41 KOPAL, Z., Numerical Analysis, London 1961. 
[l5] ABRAMOWITZ, 
M., ,”,I. 
A. (eds.), HandbookofMathematical Functions, New York 1968. 
[I61 Mom-SMITH, H. M., Phys. Fluids 3 (1960) 721. 
1171 BIRD, G. A., Molecular Gas Dynamics, Oxford 1976. 
[I81 BACHMANN, 
P., BELITZ, H. J., Report IPP 8/2 (June 1993). 
[I91 JANEV, R. K., LANGER, W. D., EVANS JR., K., POST JR., D. E., Elementary Processes in 
[20] SUCHY, 
K., in: Handbuch der Physik, Band XLIX/7, p. 57, Berlin, Heidelberg, New York, Tokyo 
[21] MITTMANN, 
H.-U., WEISE, H.-P., DING, A., ZI Naturforsch. 26a (1971) 1112. 
[22] GOLANT, 
V. E., ZHILINSKY, 
A. P., SAKHAROV, 
I. E., Fundamentals of Plasma Physics, John Wiley 
Hydrogen-Helium Plasmas, Springer-Verlag 1987. 
1984. 
Sons. New York 1980. 
