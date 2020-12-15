# GeneticAlgorithm_for_Solving_TSP
Genetic Algorithm for Solving Traveling Salesman Problem

## はじめに
[巡回セールスマン問題](https://ja.wikipedia.org/wiki/%E5%B7%A1%E5%9B%9E%E3%82%BB%E3%83%BC%E3%83%AB%E3%82%B9%E3%83%9E%E3%83%B3%E5%95%8F%E9%A1%8C) (Traveling Salesman Problem, TSP) を解くための[遺伝的アルゴリズム](https://ja.wikipedia.org/wiki/%E9%81%BA%E4%BC%9D%E7%9A%84%E3%82%A2%E3%83%AB%E3%82%B4%E3%83%AA%E3%82%BA%E3%83%A0#:~:text=%E9%81%BA%E4%BC%9D%E7%9A%84%E3%82%A2%E3%83%AB%E3%82%B4%E3%83%AA%E3%82%BA%E3%83%A0%EF%BC%88%E3%81%84%E3%81%A7%E3%82%93,%E3%83%A1%E3%82%BF%E3%83%92%E3%83%A5%E3%83%BC%E3%83%AA%E3%82%B9%E3%83%86%E3%82%A3%E3%83%83%E3%82%AF%E3%82%A2%E3%83%AB%E3%82%B4%E3%83%AA%E3%82%BA%E3%83%A0%E3%81%A7%E3%81%82%E3%82%8B%E3%80%82) (Genetic Algorithm, GA) を C 言語で実装した. 実装には [進化計算アルゴリズム入門 生物の行動科学から導く最適解](https://www.amazon.co.jp/%E9%80%B2%E5%8C%96%E8%A8%88%E7%AE%97%E3%82%A2%E3%83%AB%E3%82%B4%E3%83%AA%E3%82%BA%E3%83%A0%E5%85%A5%E9%96%80-%E7%94%9F%E7%89%A9%E3%81%AE%E8%A1%8C%E5%8B%95%E7%A7%91%E5%AD%A6%E3%81%8B%E3%82%89%E5%B0%8E%E3%81%8F%E6%9C%80%E9%81%A9%E8%A7%A3-%E5%A4%A7%E8%B0%B7%E7%B4%80%E5%AD%90-ebook/dp/B07DQFVK1H) を一部参考にした. 以下の数式を書くのに [github や note でも TeX の数式をかくぜ](https://aotamasaki.hatenablog.com/entry/2020/08/09/github%E3%82%84note%E3%81%A7%E3%82%82TeX%E3%81%AE%E6%95%B0%E5%BC%8F%E3%82%92%E6%9B%B8%E3%81%8F%E3%81%9C) を利用してみる.

## 問題の内容
n 個のノードが与えられたとき, 各ノード k (k = 1, 2, ..., n) は二次元空間にあり
![(x_k, y_k)](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%28x_k%2C+y_k%29)
はその位置を表すとする. 解は n 個のノードのうち少なくとも m 個のノードを含まなければならない. ノード k と ノード l の間の距離
![\mathrm{dist}_{(k. l)}](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cmathrm%7Bdist%7D_%7B%28k.+l%29%7D)
は  
![\mathrm{dist}_{(k, l)} = (int) (\sqrt{(x_k - x_l)^2 + (y_k - y_l)^2} + 0.5)](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cmathrm%7Bdist%7D_%7B%28k%2C+l%29%7D+%3D+%28int%29+%28%5Csqrt%7B%28x_k+-+x_l%29%5E2+%2B+%28y_k+-+y_l%29%5E2%7D+%2B+0.5%29)  
で表される. ただし ![(int) z ](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%28int%29+z+) は ![z \geq 0](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+z+%5Cgeq+0) のとき ![m \leq z](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+m+%5Cleq+z) を満たす最大の整数 m を返す. 距離は対象であることに注意.  
問題は最小距離のツアーを見つけることである. これは二次元対称 TSP と呼ばれている.

## 停止条件
アルゴリズムによって消費された CPU 時間がパラメータ timelim で与えられた制限時間を超えたときに停止する. デフォルトでは 300 に設定している.

## TSPデータについて
TSP のデータには [TSPLIB](http://elib.zib.de/pub/mp-testdata/tsp/tsplib/tsplib.html) にある .tsp ファイルを読み込む.

## コンパイル
tsp.c をコンパイルするには

```
gcc -Wall -O2 -o tsp tsp.c -lm
```

または

```
make
```

をコマンドラインから入力する.

## 実行
tsp.c を a280.tsp に対して実行するにはコンパイル後に

```
cat a.280.tsp | ./tsp
```

とコマンドラインから入力すればよい.  
実行すると以下のように出力される.

```
recomputer tour length = *****  
time for the search: ***** seconds  
time to read the instance: ***** seconds
```

1 行目は得られたツアーを評価したもので, 2 行目はアルゴリズムが消費した CPU の秒数で, 3 行目はデータファイルを読み込む時間である.  
  
また, コマンドラインから様々なパラメータを入力することもできる. 例えば, test.tsp を解く際に初期解を与えたい場合は

```
cat test.tsp test.tour | ./tsp givesol 1
```

と入力する.  
上の例では test.tour が初期ツアーのデータを書き込むファイルである. これらの初期ツアーのデータは TSPLIB 形式で書かれている必要がある. 

## tsp_view について
デフォルトではプログラムの実行後に計算されたツアーは result.tour という名前のファイルに出力される. result.tour を tsp_view ディレクトリに移動し  

```
cat result.tour | ./tsp_view.py
```

とコマンドラインから入力するとツアーの画像ファイルが出力することができる. (gnuplot のインストールが必要)

## cpu_time.c について
プログラム cpu_time.c は土村展之先生が作ったものであり, これを使い計算時間を推定している. 最新版は [土村先生のホームページ](http://tutimura.ath.cx/~nob/c/) で公開されている.

## 工夫した点
### 順序配列
GA では  order _ representation 関数で順序配列を用いて遺伝子の操作を行った.

### 二点交叉
tow _ point _ crossover 関数で二点交叉を実装した.
