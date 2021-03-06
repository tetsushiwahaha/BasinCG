#basinCG

<img align="center" src="https://user-images.githubusercontent.com/52724526/85823990-b020e500-b7b9-11ea-8064-abbd4f4615ae.png">

「CによるカオスCG」オマケ <棲み分け図> 作成のしおり 

* 1st.	Sat Jun 25 15:00:01 JST 1994
* rev.	Thu Sep  8 13:30:27 JST 1994
* rev.	Wed Jan  3 17:17:22 JST 1996
* rev.	2000年 7月14日 金曜日 10時16分41秒 JST
* rev.	2010年 10月15日 金曜日 19時38分21秒 JST
* rev. 2020/06/26

1. 内容

このディレクトリには以下のファイルがあります．

* `README.md`: このファイルです．
* `Makefile`: `Basin`, `BasinPutColor` を作成するためのmakefileです．
* `BasinDefs.h`: `basin.c`で用いるデータを定義しています．
* `Basin.c`: 第一引数で与えられたファイルの記述にしたがって，棲み分けを
計算するプログラムです．このプログラムでは，
`CによるカオスCG'のプログラム例4.1の漸化式を計算させています．
結果として，第二引数で示されるファイルに
ヘッダ，周期点，収束回数データを書き出します．
* `BasinPutColor.c`:    第一引数にbasinプログラムによって得られた
出力ファイルを指定します．そのデータを対話形式で色づけ
できます．
* `in0`, `script`: サンプル設定ファイルです．4章を参照してください．


2. プログラムの作り方

画像フォーマットとしてPNG を用いますので，関連するライブラリをインストー
ルする必要があります．MacOS ならば，macports を導入していれば，

    % sudo port install libpng 

だけでOKです．macportsは /opt/local 以下にライブラリなどもインストールさ
れますので，-lz もそこを指すように Makefile では書いてます．
そのあと

    % make all

とします．コンパイル時には，いくつかの
警告が出ますが，無視して構いません．
これで `Basin` および `BasinPutColor`が得られます．
UNIXプラットフォームによっては若干`Makefile`を手直しする必要が出て来るかも
知れません．詳細は `Makefile` を確認ください．


3. 使用方法

3.1 引力圏計算プログラム `Basin`

3.1.1 設定ファイルの記述

 設定ファイルでは，パラメータの値，x-y平面の大きさ，メッシュの粗さ(ピクセル)
 写像回数，発散とみなす条件，精度1, 2 を記述します．

```
(例)
0.6 0.98 4.0			/* a, b and c */
-40.0 40.0			/* xmin xmax */
-40.0 40.0			/* ymin ymax */
50 50 				/* width and height */
2000				/* mapmax */
1000.0				/* divergence */
1e-6				/* eps */
1e-5				/* eps */
```

 C言語ふうにコメントが書き込めます．この例では，

| 行 | 意味 | 例の値|
|----|------|-------|
| 1  |パラメータ| `a = 0.6, b = 0.98, c = 4.0` |
|2 | x-y平面             | (-40,-40) 〜 (40, 40) |
|3 | ウィンドウの大きさ  | 50x50 |
|4 | 最大写像回数        | 2000回 |
|5 | 発散とみなすノルム  |  1000.0 |
|6 | 周期点の判定精度    | 1e-6|
|7 | 棲み分けの判定精度  | 1e-5|

という意味になっています．これらの記述順番は変えられません．
計算時間のポイントとしては，

1. パラメータにより計算時間が大幅に変化する．
1. 写像回数の多い周期点が多いほど計算時間がかかるこの例の場合，
bが1.0に近いとカオティックな応答になるため，写像回数が
増え，計算時間がかかる．
1. ウィンドウの大きさに比例(非線形比例)して時間がかかる．

が挙げられます．最大写像回数を越えても収束しない解はカオスとみなすと
よいでしょう．
なお，棲み分けの判定精度を周期点の判定精度より同じか，厳しく(小さく)
すると正しい棲み分け図が得られません．

3.1.2 使用方法

    % ./Basin in 

inという設定ファイルを読み込み，出力として，in.baz ファイルを生成します．
baz ファイルは引力圏のデータを zlib を用いて圧縮してあるファイルです．
サンプルに，in* の名前でいくつか設定ファイルを置いてあります．


3.2 画像加工プログラム BasinPutColor

BasinPutColor は Basin プログラムで生成されたファイルを解析して，対話的に
色情報を入力して，png ファイルに出力するプログラムです．

3.2.1 設定ファイル

設定ファイルはありません．できたら，設定ファイルを使用するように
作り変えていただきたいぐらいです．

3.2.2 使用方法

    % ./BasinPutColor in 

とすると，カレントディレクトリの in.baz というファイルを
探し，存在するなら内容を解析し，以下のように表示します．
```
header length = 40
type of period = 5
ix = 50
iy = 50
id:0 period:6 freq:1532 
id:1 period:1 freq:379
id:2 period:8 freq:220
id:3 period:1 freq:368
id:4 period:1 freq:1
```

ここで，これらの意味は，

| 記号 | 意味 | 
|------|------|
|header length  | これは気にしなくて構いません．|
|type of period | 周期点の種類 |
|ix, iy         | ウィンドウの大きさ |
|id             | 周期点の番号．発見した順にid番号が付いている． |

などです．この例では，

周期点が5種類で，それぞれ，
```
6周期点が1532点
1周期点が379点
8周期点が220点
1周期点が368点
1周期点が1点(ただし，これはサドルである)
```
という意味です．全部の点を加算すると 2500点(= 50x50)になっているはずです．

続いて，

    beta =

と聞かれます．これはべタ塗りするかどうかの指定で，0なら
濃淡塗りをし，1ならベタ(単一色で棲み分け)となります．
これ以降も yesなら1, noなら0を入力します．

    automatic color assigning =

自動色塗りを行なうかどうか，です．1を入力すると
勝手な色が割り付けられます． (だったと思う)
滅多に使わなくなりましたから，もう忘れました．(^_-)

    using another color for every period = 

周期点に出会った時，それぞれの周期点に属する引力圏を
別の色で塗り分けるかどうかを聞かれます．1を入力した
場合，例えば6周期点ならば，このアトラクタに引き込まれる
領域を6色の領域で塗り分けます．
つまり，一旦6周期点に収束したと判定された初期値について，
もう6周期に落ち着くという情報を得たわけですから，再び同じ
初期値から，今度は6回おきに写像の値を取っていくと，
6周期点のbasin of attractionsが取れるわけです．
この後にその色を指定せねばなりません．

    add automatic bias =

収束回数でピクセルの発色に重みをつけるのですが，
デフォルトでは全ての点の収束回数を一度スキャンし，
一番収束回数の少ない(もしくは多い)点を最低輝度に
するようになっています．(値0〜1)．輝度を少なく
調整したい場合に用いたと思います．
1を入力すると，

    bias =

と聞かれますので，値を指定します．
値は0〜1の間の小数です．単純にその値が
データ値(0〜1)と積がとられ，新たなデータ値となります．
次に

    addbias =

となって， automatic biasを入力しても，またバイアスを加えるかどうか
聞かれてしまいます．(トホホ) 1を入力すると，

```
period:6
min[0] = 1295
max[0] = 1721
bias =
```

とその周期点の最大最小値を表示し(データとしては最小値が0で最大値が
1になるように正規化されている)，それらにかけるバイアス値を指定します．
ここまで書いて，こういった使い方からいうとバイアスという言い方は
適当でないという気がしてきました．小数0〜1の値をかけるということは
全ての値は小さくなってしまいますね．たしか，この値を後に２乗して
使ったりしましたが，忘れました．どなたか，ここらへんをスッキリ
整理してバージョンアップを図ってくれませんか．

まぁ automatic bais, addbiasともに0を答えておけば深みにハマりません．

    period 6: id:0-0, f(1532): h s i=

automatic colorを1としない限り，このように周期点に割り当てる
色を聞いてきます．

色は色相(hue)，飽和度(sturation)，輝度(intensity)で構成されます．

色相は角度になっており，
```
0 : 赤
60 : 黄色
120 : 緑
180 : シアン(水色)
240 : 青
300 : マゼンタ(赤紫)
```
っと360度回っています．

飽和度は1: 原色，0: 白 の間の値をとります．

輝度は，1: 明るい，0: 暗い(黒)の間の値をとります．

それで，
`period 6: id:0-0, f(1532): h s i=`
の`h'は0〜360の角度，`s'は0〜1の実数を与えます．これは
色相と飽和度です．
最後の`i'は収束回数の大小に対して，明るさの大小をどう割り振るか
を決めます．0ならば，収束回数の最も少ない点は一番暗く，
1ならば収束回数の最も少ないものを一番明るく塗ります．

たとえば，

    period 6: id:0-0, f(1532): h s i= 240.0 0.5 1

は6周期点のある引力圏を，青のパステル(中間色)で塗りますが，
収束回数の一番少ないもの(すなわち周期点そのもの)を一番
明るくなるように塗ります．

    -1 : hue & sat & in =

カオスとみなされた点，発散とみなされた点，準周期とみなされた点など
は収束回数を-1, -2, -3に割り当てて記憶しています．それらの色を
どう与えるかです．

    -1 : hue & sat & in = 0 0 0

とすると，その部分は真っ黒になります．
Basin.cでは-1はカオスを表します．
他の棲み分けプログラムでは-2 や -3 を返すものがあり，
カオス，準周期，発散なんかを表したような気がします．

以上がBasinPutColorに与えるデータです．計算結果はただちに png フォーマッ
トのファイルとして保存されます．



4. サンプルファイル

さて，このディレクトリには in*, script*というサンプルも
入っています．

●サンプルオペレーション

    % ./Basin in0 

と入力するとin0に従って引力圏の計算を開始します．
続いて，

    % ./BasinPutColor in0 < script0

とするとデータout0と共にout0.rasが得られます．script0は
実際は対話的に与えるデータを面倒なので，手入力する内容を
ファイルにしたまでの話です．

in0.png の表示には例えば

    % display in0.png

としてください．

もう一つ，in3というファイルも用意していますが，カオティックな
過渡応答が長いパラメータなので，相対的に計算に時間がかかります．

5. バグ

バグだらけです．なんのメンテもしていません．
なんとかしてくれません？
