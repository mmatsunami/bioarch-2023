# フェージングハンズオン

河合　洋介（国立国際医療研究センター）

2023年12月06日作成

## 準備

### サーバへのログイン

まずsshコマンドを用いてサーバにログインしてください。

"user"と"node_name"は割り当てられたものを指定して下さい。

```sh
ssh user@gwa.ddbj.nig.ac.jp
ssh node_name
```

設定ファイルをよみこみます。

```
source /home/bioarchaeology-pg/setup.sh
# Javaが実行できるか確認
java --version
```

### 使用するプログラム

[BEAGLE5.4](http://faculty.washington.edu/browning/beagle/beagle.html)

今回はスパコンにダウンロード済みのファイルを使います。

#### プログラムの確認
```
java -jar /home/bioarchaeology-pg/kawai/beagle.22Jul22.46e.jar
```

ヘルプメッセージが表示されることを確認

### 使用するデータ
シミュレーションで作った250個体の*Fictus yaponesiae*のゲノムデータ(`yaponeasi_SP_250.vcf.gz`)を使います。

#### データの確認
```
# ファイルの存在の確認
ls -l /home/bioarchaeology-pg/kawai/yaponesia_SP_250.vcf.gz
# ファイルの中身の確認
bcftools view /home/bioarchaeology-pg/kawai/yaponesia_SP_250.vcf.gz | less -S
```
ここでVCFファイルのフォーマットを確認する。特に遺伝子型が`0/1`のように`/`で区切られていることに注目する。

### 解析ディレクトリの作成
```
# ホームディレクトに移動して作業ディレクトリを作成
cd
mkdir phasing
cd phasing
```

## 実行
### BEAGLEでフェージングを実行
```
java -jar /home/bioarchaeology-pg/kawai/beagle.22Jul22.46e.jar gt=/home/bioarchaeology-pg/kawai/yaponesia_SP_250.vcf.gz out=yaponesia_SP_250.phased nthreads=8
```
この例では`yaponesia_SP_250.vcf.gz`を入力ファイル(`GT=`)にして、`yaponesia_SP_250.phased`をプレフィックスに出力する(`out=`)。`nthreads=`で実行に使用するCPUの数を指定する。
終了後にフェージング結果（`yaponesia_SP_250.phased.vcf.gz`）がカレントディレクトリにできていることを確認する。

### フェージング結果の確認
```
bcftools view yaponesia_SP_250.phased.vcf.gz | less -S
```

先ほどの入力データとの違いを確認する。

### ヒトゲノムデータの確認
国際1000人ゲノム計画の日本人集団(JPT)の全ゲノムシークエンス解析で得られたフェージング済みのVCFを確認します。
`/home/bioarchaeology-pg/kawai/JPT`

```
bcftools view /home/bioarchaeology-pg/kawai/JPT/JPT.phased.chr1.vcf.gz | less -S
```
後でIBDNeを使った人口動態推定に使います。




