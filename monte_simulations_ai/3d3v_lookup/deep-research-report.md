# D + D+ 中性粒子 Monte Carlo における elastic energy transfer 不一致の文献ベース原因解析

## 背景と観測された不一致の位置づけ

ご提示の状況は「**同じ（または同等と主張している）弾性散乱データ**」から、  
(A) 角度サンプリングによる **lab-frame の逐次 ΔE（collision / track-length）集計**と、  
(B) **EIRENE/DEGAS 系の transport-moment テーブル（I\_1\_0, I\_1\_1, I\_1\_2）による track-length 推定**  
が、**CX では整合するのに EL（elastic）で大きく乖離する**という問題です。

文献側（EIRENE/DEGAS 系）では、**弾性衝突の運動量・エネルギー交換は「微視的散乱角分布」から導かれる“輸送（transport）断面”のモーメント**として整理され、Monte Carlo ではそれを **collision estimator（CL）** で直接サンプリングしても、**track-length estimator（TL）** でテーブル参照しても、**同一の衝突積分（collision operator）の平均**を狙う、という構図です。実際、EIRENE 系の文献では **CL と TL の整合性テストを“実装整合”の重要なチェックとして位置づけ**ています。citeturn41view5

したがって、まず押さえるべきは不一致が

- **（i）本質的に「同じ量」を比べていない**（frame・符号・定義のすれ違い）
- **（ii）同じ量を比べているが、どちらか（特にテーブル生成側）が物理・数値的に破綻**
- **（iii）同じ量で両者とも正しいが、D + D+（等質量・ほぼ熱平衡）ゆえに “net power が小さく相対差が増幅”**

のどこに属するか、です。以下、文献で定義されている「I テーブルの意味」と「一致条件」「代表的な落とし穴」を整理し、そのうえで **apples-to-apples の検証手順**を提示します。

## Transport-moment と “I\_1\_0 / I\_1\_1 / I\_1\_2” が表すもの

### 導出の核は「輸送（momentum-transfer）断面 σ^(1)」である

EIRENE 系の理論整理では、弾性衝突における運動量・エネルギー交換は、単なる全断面積（total cross section）ではなく、**散乱角で重み付けした輸送断面**（典型的には momentum-transfer / diffusion cross section）で決まります。

entity["people","Vladislav Kotova","solps-eirene author"]（B2-EIRENE のアップグレードに関する報告書）では、いわゆる “Maxwellian molecules” の議論の中で、**σ^(1)** を次の形で提示しています（χ は **CM 系の deflection angle**、v\_r は相対速度、b は impact parameter）：  
**σ^(1)(v\_r) v\_r = 2π v\_r ∫ [1 − cos χ(v\_r, b)] b db**。citeturn41view0  

この式は、「散乱角分布（微分断面）を（1−cosχ）で重み付けして積分した量が、運動量拡散（=輸送）を支配する」ことを明示します。従って、**total の反応率（reaction rate）が合っても、σ^(1)（=輸送断面）が合うとは限らない**、という構造が文献の定義から直ちに出てきます。citeturn41view0

### I(1,0), I(1,2)…は「σ^(1) を Maxwell 背景で畳み込んだモーメント（速度・エネルギー重み付き積分）」

Kotov の同報告書では、σ^(1) v\_r が定数なら衝突積分が簡単になり、そこから運動量交換・エネルギー交換の閉形式を導ける一方、現実の σ^(1)(v\_r) は速度（エネルギー）依存なので、**その速度依存を保持した積分（文中の “integrals”）を扱う必要がある**と述べています。citeturn41view4turn41view5

そして「定数（K\_m）の代わりに I^(1,0)(E,T) を使う」として、I^(1,0) と I^(1,2) が、運動量・エネルギー交換の式に入ることを示しています。citeturn41view4turn41view5  
この I^(1,0), I^(1,2) はまさに、**σ^(1)（角度重み付き）を背景 Maxwell（温度 T）で畳み込んだ“輸送モーメント”**です（E はテスト粒子のエネルギーをどの frame で評価するかが重要になります。後述）。citeturn41view4turn41view5

ご質問の I\_1\_0, I\_1\_1, I\_1\_2 は記法差はあり得ますが、EIRENE/DEGAS 系で一般に  
**第1添字 “1” が σ^(1)（momentum-transfer / diffusion）に対応する次数**、  
**第2添字が速度（エネルギー）重みのモーメント次数**  
を表す系列で、Kotov の I^(1,0), I^(1,2) と同系列の量として理解するのが自然です。citeturn41view0turn41view4turn41view5

### D + D+（等質量・背景流速0）で重要な“キャンセル構造”が、文献の式レベルで確認できる

Maxwellian molecules の簡約式では、全体のエネルギー交換率が **〈E₂〉−〈E₁〉**（と drift 項）に比例する形が露出します：  
(dE/dt)\_{21} ∝ (〈E₂〉 − 〈E₁〉 + (m₁−m₂)/2 (u₁·u₂))。citeturn41view1  

ここで m₁=m₂、背景 drift が無い（u₁·u₂ が消える/小さい）なら、**“エネルギー差だけ”が駆動**になります。  
したがって D + D+ のような等質量で、かつテスト中性が背景イオンと“ほぼ同じ温度（2 eV 近傍）”に寄っている区間では、**net の EL power が本質的に小さくなり、推定量の微小なバイアスが相対的に巨大化**し得ます。citeturn41view1  

この “小さい差分を見ている” 構造は、後述する「full simulation で 10–100 倍に増幅される」観測と整合的です。

## Naive 角度サンプリングと I テーブルは理論的に同じ量か

### 結論：狙っているのは同じ “衝突積分の平均” だが、一致には明確な条件がある

文献側の立場では、CL（個々の衝突をサンプリング）と TL（経路に沿って平均寄与を積算）は、同じ衝突モデルを実装する限り **一致すべき**であり、整合テストは「衝突実装の一貫性チェック」とされます。citeturn41view5  
つまり理論的には、**I テーブルによる lookup energy transfer は、角度サンプリングで得る lab-frame ΔE の期待値（適切な重み付き平均）と同一**を狙います。

ただし、一致のために最低限必要な条件が、EIRENE の公式マニュアルにも明記されています。EIRENE では背景衝突相手の速度は単なる Maxwell サンプルではなく、**（必要に応じ drift を含む）cross-section weighted Maxwellian**からサンプルし、断面は **相対速度で評価**すると説明されています。citeturn41view11turn41view13  
この “σ(v\_rel)·v\_rel で重み付けされた相手サンプル” が崩れると、**reaction rate だけでなくエネルギー交換も意味が変わる**ため、apples-to-apples ではなくなります。citeturn41view13

### 一致条件を「検証可能な形」に落とすと

角度サンプリング法（A）が lookup（B）と一致するための、文献に整合するチェック観点は次の通りです：

- **同じ衝突演算子（同じ differential cross section / transport cross section）を使っていること**  
  I テーブルは σ^(1)（= 角度重み付き輸送断面）を核に構築されるため、A 側が用いている角度 CDF が **同じ χ（CM deflection angle）定義**と同じ dσ/dΩ 正規化に基づく必要があります。citeturn41view0  
- **背景イオン速度サンプリングが “σ(v\_rel)·v\_rel 重み”を満たすこと**  
  EIRENE マニュアルは、衝突相手速度分布が σ(v\_rel)·v\_rel により選択されること、またそれを rejection で実現していることを述べています。citeturn41view13  
- **同じ frame（R-frame / L-frame）で E を定義し、同じ量（テスト粒子の ΔE か、背景の ΔE か、符号規約）を比較すること**  
  Kotov は背景ドリフトがある場合に、R-frame（背景流速0の frame）で計算した量を L-frame に変換する手順を明示しています。citeturn41view4  
  ここがずれると、A の “lab-frame ΔE” と B の “（R-frame で定義された）エネルギー源” を比較してしまい、見かけ上大きな差が出ます。citeturn41view4  

以上を満たすなら、**「I モーメントでの lookup energy transfer」と「角度サンプリングでの平均 ΔE」は同じ量**とみなしてよい、というのが文献整合の答えになります。citeturn41view5turn41view11turn41view13

## 同じ Bachmann 系データでも「CDF と transport-moment が大きくズレる」代表原因

ここはご観測（reaction_rate は合うが I テーブルが大きくズレる）に直結します。文献から引ける“ズレやすい本質”は大きく3系統あります。

### 角度分布の僅かな不整合が、輸送断面（σ^(1)）を大きく誤らせる

輸送断面 σ^(1) は (1−cosχ) の重みを持つため、**全断面積（反応率）とは角度領域の寄与が違う**、ということが決定的です。citeturn41view0  

entity["people","Predrag S. Krstic","ornl physicist"] と entity["people","David R. Schultz","ornl physicist"] の短報は、輸送断面の数値積分が難しい理由をかなり明確に述べています：

- momentum transfer（拡散）断面は **CM 後方角**の寄与が支配的、  
- 粘性断面は **中間角**が支配的、  
- したがって dσ/dΩ を **0–180° 全域で高精度に保つ必要がある**が、前方と後方の dσ/dΩ が **10^5** も異なる場合があり、これが数値的難しさの核心、  
としています。citeturn35view0  

この構造の帰結として、

- 全断面（→反応率）は合っているのに
- σ^(1) やそれを畳み込んだ I モーメントが大きくズレる

という現象は、**物理的にも数値的にも起こり得る**、というのが文献からの含意です。citeturn35view0turn41view0  
特に **forward peak が強い**データでは、全断面は forward で稼げても、(1−cosχ) 重みで forward が抑制され、**“尾部（大角度側）”の誤差が輸送断面を支配**する、という典型的な増幅経路が成立します。citeturn35view0turn41view0  

### “cumulative scattering angle” の角度定義・測度が一致していない

輸送断面の定義に現れる角度 χ は **CM 系 deflection angle** です。citeturn41view0  
一方、CDF ファイルが持ちうる角度表現には少なくとも次があり、取り違えが起きると輸送モーメントが破綻します：

- 角度が χ（CM deflection）なのか、別の “scattering angle” なのか（用語揺れ）
- CDF が θ に対する CDF なのか、cosθ に対する CDF なのか（測度が違う）
- “dσ/dΩ” を基にした分布なのか、“dσ/dθ” を基にした分布なのか（sinθ 因子の有無）

文献側は、少なくとも σ^(1) に関して「(1−cosχ) を用いる」ことと「χ が CM deflection」であることは明確です。citeturn41view0  
したがって、**CDF 側の角度がこの χ と同じでない**場合、A と B の一致は原理的に壊れます。

ご観測の「reaction_rate は整合するが I テーブルがズレる」は、角度分布の取り違え（あるいは sinθ 因子欠落）と非常に相性が良い症状です。反応率は角度に依らないので破綻が露見しにくく、輸送モーメントだけが壊れます。citeturn41view0turn35view0  

### エネルギー変数 E の取り扱いが一致していない（E_r vs E_lab）

Kotov は弾性衝突データのフィッティングで、相対運動エネルギー E_r を、便宜上 moving frame の E_lab に置き換えていること、そして両者の関係を式として明示しています。citeturn41view2  
さらに、量子計算との比較で「古い比較が低エネルギーでも大きくズレたのは **E_lab → E_r の rescale をしていなかったため**」という指摘をしています。citeturn41view2  

つまり、Bachmann 系データの利用では **“どのエネルギーを独立変数として tabulate / fit しているか”** が一致条件の一丁目一番地で、ここを誤ると I テーブル側の再積分や lookup が大きく崩れ得ます。citeturn41view2turn41view4  

## 文献に現れている注意点・不一致報告との対応付け

### Bachmann 1995（Bachmann–Reiter）に起因する注意点

Kotov は、参照文献 [44]（文献リスト上、entity["people","P. Bachmann","fusion edge atomic data"] と entity["people","Detlev Reiter","fusion plasma code developer"] の 1995 年論文に対応）の式について、**RE（エネルギー交換率）の式に誤りがある**と明記しています。citeturn29view4  
（このページはスクリーンショット取得が一時 503 で落ちたためテキストビュー引用ですが、同報告書の該当箇所の記述です。）

この種の「元論文の式の誤植・取り違え」は、そのままデータ生成（I テーブル生成）コードへ伝播すると、**反応率（全断面）は合うが、エネルギー交換だけが壊れる**という形のバグになり得ます。

### EIRENE 側が公式に明示している “衝突相手サンプルの重み” と “CX のデフォルト仮定”

EIRENE マニュアルは、衝突に入る“二次粒子”（あるいは衝突相手）の速度ベクトルを **cross-section weighted Maxwellian**からサンプルし、断面は相対速度で評価する、と書いています。citeturn41view11turn41view13  
ここは naive 実装で最も落としやすい（が理論的に “同じ量” を保証する要件）です。

また、共鳴 CX のデフォルトモデルについて、

- CM 系散乱角 π
- exothermicity 0（CM 系で運動エネルギー保存）
- “exchange of identity” として速度を保ったまま電荷を入れ替える（という運動学モデル）

を明示しています。citeturn41view12  

このモデルは、まさに CX で naive と lookup が整合している、というご観測と整合的です（CX は“モデルが強く規定されている”ため、実装揺れが小さい）。citeturn41view12turn41view5  

### “EL はズレやすい” こと自体は、B2-EIRENE 文献でも示唆されている

Kotov は CL と TL の比較を行い、**CX では簡略化しても差が 25% 以内**だが、**弾性（例：分子–イオン弾性）では factor 2 以上ズレ得る**と述べています。citeturn41view5  
これはそのまま「弾性の energy transfer は実装・近似・テーブル生成に敏感」という含意で、CX だけが合う/EL だけがズレるという今回の現象の“方向性”に一致します。citeturn41view5  

ただし、あなたのケースは factor 10–100 とさらに大きいので、**単なる近似誤差ではなく、定義・frame・モーメント生成のどこかが破綻している可能性**が高い、と推定するのが自然です（根拠は次節の apples-to-apples 検証で切り分け可能です）。citeturn41view5turn41view0turn35view0  

## Apples-to-apples の比較設計

ここでは「どの量同士を比べればよいか」を、文献の定義に合わせて“手順化”します。ポイントは、反応率（σ^0）と輸送モーメント（σ^(1) とその速度畳み込み）は別物なので、**両方を同じ入力データから再構成できるか**をテストすることです。citeturn41view0turn35view0  

### 比較対象の明確化

まず、比較したい量を次のどれにするか固定します（ここがズレると永久に一致しません）：

- **(i) テスト中性粒子の平均エネルギー変化率**（d〈E\_t〉/dt あるいは d〈E\_t〉/ds）
- **(ii) 背景イオン（bulk）の平均エネルギー変化率**
- **(iii) コードが出力する “EL power”（符号規約つき）**

EIRENE マニュアルは “bulk ion energy” と “collision kinetics” の区別・オプションの枠組みを示しており、同じ衝突でも「どの粒子群のエネルギーとして計上するか」が設定依存になり得ることが読み取れます。citeturn41view10turn41view11  
したがって、まずは **lookup 側の dtr\_eng が “誰のエネルギー損得” を表現しているか**（テスト粒子か bulk か）を確認し、それに合わせて naive の集計量（−ΔE の符号を含む）を定義し直してください。

### 単セル（0D）での3経路照合が最短

次に、空間輸送を絡めない 0D/単セルで、同一条件下で3つの経路を照合すると原因が局所化します。

- **経路1：数値積分（deterministic quadrature）**  
  背景イオン速度分布は、EIRENE が要求するように σ(v\_rel)·v\_rel 重み付きである必要があります。citeturn41view13  
  角度は CDF に従って（ただし χ が CM deflection であること）。citeturn41view0  
  このとき、期待値としての 〈ΔE〉（あるいはエネルギー交換率）を高精度に算出します。  
- **経路2：naive Monte Carlo（あなたの方法A）**  
  EIRENE の記述どおり、衝突相手速度分布が σ(v\_rel)·v\_rel で重み付けされていることを再現します（rejection など）。citeturn41view13  
  そのうえで lab-frame ΔE を統計平均し、経路1と一致を確認します。  
- **経路3：I テーブル lookup（あなたの方法B）**  
  同じ E（どの frame かを固定）と Ti で I\_1\_0, I\_1\_1, I\_1\_2 を引き、同じ量（テスト or bulk）に変換して比較します。Kotov は R-frame → L-frame の変換構造を明示しているので、**lookup 側がどの frame で定義されているか**をここで確定できます。citeturn41view4turn41view2  

この3経路のうち、

- (1) と (2) が一致するが (3) がズレる → **I テーブル生成か、I の定義/単位/適用式が違う**
- (2) が (1) とズレる → **角度 CDF の測度（sinθ）や χ 定義、相手サンプルの重み付けが違う**
- (1) と (3) は一致するが full simulation の time-average でだけ爆発 → **“net power が小さい”ことによる増幅（次節）＋空間輸送との相互作用**

に分類できます。

### “reaction_rate は合うが I が合わない” を再現するための最小チェック

ご観測（反応率は整合、I が不整合）をそのまま診断変数にすると、最小チェックは次です。

1) **全断面 σ^0(E) の同一性**：これは既に合っているとのことなので OK。  
2) **輸送断面 σ^(1)(E) の再構成**：  
　CDF（χ）から 〈1−cosχ〉 を計算し、σ^(1)=σ^0·〈1−cosχ〉（※角度定義が一致している場合）を検査。ここがズレるなら、CDF と σ^(1) 定義が違います。σ^(1) の定義自体は文献で明確です。citeturn41view0  
3) **E の変数変換**：E\_lab と E\_r の取り違えがあると低エネルギーほど致命的、という指摘が文献にあります。citeturn41view2  
4) **Bachmann 1995 の式誤りの影響**：RE の式誤りが指摘されているため、I テーブル生成コードがその式に依存していないかを確認。citeturn29view4  

この 2) が最重要です。反応率は 0 次モーメント（角度無重み）ですが、I は事実上 **角度重み付きモーメント**なので、そこでズレるのは “典型症状” です。citeturn41view0turn35view0  

## Naive TL を真値とみなしてよいか

### 「理論真値」になり得る条件と、危険な点

naive 法（角度サンプリングで lab-frame ΔE を積算）は、実装が以下を満たすなら **衝突積分の直接 Monte Carlo**として理論的に正しく、テーブル法の参照（ground truth）になり得ます。

- 衝突相手速度が **σ(v\_rel)·v\_rel 重み付き**になっている（EIRENE はこれを明示し、rejection で実現）。citeturn41view13turn41view11  
- 角度分布が σ^(1) の定義と整合する χ（CM deflection）に基づく（σ^(1) の定義は (1−cosχ)）。citeturn41view0  
- 比較する “エネルギー” が同じ frame / 同じ粒子群 / 同じ符号規約で揃っている（R-frame → L-frame の注意）。citeturn41view4turn41view2  

一方で、naive 法には今回の条件（D + D+、Ti=2 eV、背景流速0）で特有の危険が2つあります。

1つ目は **差分のキャンセル**です。等質量・低ドリフトでは、エネルギー交換は 〈E₂〉−〈E₁〉に比例する形で小さくなり得ます。citeturn41view1  
このとき、naive の “衝突ごとの ΔE” は分散が大きいのに平均が小さいため、統計誤差・わずかなバイアスが **net power を 10–100 倍に見せる**ことが現実に起きます（あなたの観測と整合）。citeturn41view1  

2つ目は **輸送断面の高感度**です。dσ/dΩ の角度依存が極端な系では、輸送断面（そして I モーメント）が「特定角度帯の寄与」に強く依存し、角度 CDF の表現誤差が増幅されます。citeturn35view0turn41view0  
このケースでは、naive のほうが“真値”というより「入力 CDF の誤りを忠実に再現する」だけになり得ます。

従って結論としては：

- **単セル（0D）で deterministic quadrature と一致することを確認できた場合に限り**、naive TL を参照値とみなすのが妥当  
- full simulation の time-average だけでの 10–100 倍差は、まず「キャンセル増幅」か「frame/符号不一致」を疑い、次に「I テーブル生成の破綻」を疑う

が、文献整合の評価になります。citeturn41view1turn41view4turn41view0  

---

### 追加で強く疑うべき点（ご観測に特化）

最後に、ご観測のうち「dd\_00\_elastic.cdf は reaction\_rate は整合するが、I\_1\_0/I\_1\_1/I\_1\_2 が同じ sigma+angle から再積分した値とかなりずれる」という症状に対して、物理というより “定義・生成” 観点での優先疑い順をまとめます（どれも上の apples-to-apples テストで判定可能です）。

- **σ^(1) を作るときの角度測度（sinχ）取り違え**：反応率は救われるが輸送モーメントが壊れます。σ^(1) は (1−cosχ) で重み付ける定義。citeturn41view0  
- **χ（CM deflection）ではない角度の CDF を χ として扱っている**：輸送断面が一致不能になります。citeturn41view0  
- **E\_lab / E\_r の取り違え**：低エネルギーで壊滅的になり得ることが文献で明示されています。citeturn41view2  
- **Bachmann 1995 の RE 式誤りの伝播**：I テーブル生成にその式を使っていると “反応率は合うがエネルギー交換が壊れる” 形になります。citeturn29view4  
- **比較している量が違う（bulk vs test、R-frame vs L-frame、符号規約）**：Kotov は frame 変換の式構造を明示しており、ここが揃っていないと見かけ差が出ます。citeturn41view4turn41view1  

以上を踏まえると、現状の情報だけからのもっとも“文献に沿った”推論は、

- **(A) と (B) は原理的に同じ量を狙う設計だが**、  
- D + D+ の低エネルギー弾性は **輸送断面の高感度**と **等質量・熱平衡近傍のキャンセル**のため、  
  “反応率は合うのにエネルギー交換だけが壊れる” という症状が起きやすく、  
- さらに Kotov が指摘する **Bachmann 1995 のエネルギー交換式の誤り**や、E\_lab/E\_r の rescale の取り扱いは、I テーブル生成コード側の系統誤差要因として強く疑うべき、  

という整理になります。citeturn41view0turn35view0turn41view1turn41view2turn29view4

それはかなり重要な切り分けです。

**Collision estimator と Track-length estimator で EL が一致している**なら、少なくとも **Monte Carlo の estimator 形式そのもの**（CL か TL か）が主因である可能性はかなり下がります。EIRENE 系でも、CL と TL は同じ衝突演算子の平均を評価する実装として扱われ、整合は実装健全性チェックの一つです。

この場合、疑う順番はこうなります。

1. **比較している物理量が本当に同じか**
   test neutral の ΔE なのか、bulk ion への power なのか、符号規約は同じか、R-frame/L-frame の混同はないか、です。EIRENE/DEGAS 系の transport-moment 式は frame と定義量に敏感です。

2. **入力の角度分布と transport moment が同じ collision kernel を表しているか**
   反応率が合っても、transport cross section は別です。EIRENE 系の弾性輸送は角度重み付きの
   [
   \sigma^{(1)} \propto \int (1-\cos\chi), d\sigma
   ]
   に支配されるので、全断面積が合っていても、角度分布の定義や測度が少し違うだけで (I_{1,0}, I_{1,1}, I_{1,2}) は大きくずれ得ます。

3. **table 生成側の問題**
   あなたの観測どおり、**reaction rate は整合するのに I tables が再積分値と合わない**なら、最優先で疑うべきは estimator ではなく **I-table 生成法**です。特に

   * scattering angle の定義（CM deflection angle か）
   * CDF が (\theta) か (\cos\theta) か
   * (d\sigma/d\Omega) と (d\sigma/d\theta) の取り違え
   * (E_{\rm lab}) と relative energy (E_r) の取り違え
     が典型です。Kotov の報告では、低エネルギー比較で (E_{\rm lab}\to E_r) の rescale 不備が問題になることも明示されています。

4. **D + D(^+)・等質量・低温での差分増幅**
   等質量・ほぼ熱平衡では net EL power は差分量になりやすく、絶対誤差が小さくても相対差は大きく見えます。これは依然として妥当です。なので full simulation で 10–100 倍差が見えても、まずは「真の巨視的加熱差」ではなく、**小さい差分への系統誤差の増幅**を疑うのが自然です。

なので、いまの情報を踏まえると結論はこうです。

**結論**

* **EL の CL と TL が一致している** → Monte Carlo estimator のせいではなさそう。
* 真の争点は **「naive angle-sampling kernel」と「I-table kernel」が本当に同じか」** です。
* とくに **dd_00_elastic.cdf で reaction_rate は合うが I tables が合わない** なら、物理的には **transport cross section / transport moment の生成・定義ミス** を第一候補にすべきです。

apples-to-apples の最短チェックはこれです。

* 同じ (E), (T_i) で

  1. CDF から (\langle 1-\cos\chi\rangle) を直接出す
  2. (\sigma^{(1)}=\sigma^0\langle 1-\cos\chi\rangle) を作る
  3. それを Maxwell 背景で畳み込んで (I_{1,0}, I_{1,1}, I_{1,2}) 相当量を再構成する
  4. ファイル内 I table と比較する
* ここで崩れたら、**Monte Carlo ではなく table / CDF / frame 定義**の問題です。

この新事実を入れると、私は **「naive TL が真値か」よりも「naive kernel と lookup kernel が同一か」** を中心に見るべきだと評価を更新します。
必要なら次に、**あなたの式の (sp0, dtr_eng) が test-particle energy と bulk-ion power のどちらを返しているか**を、文献式に沿って1項ずつ対応付けて整理します。
