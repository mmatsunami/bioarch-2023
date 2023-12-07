# 集団サイズ推定ハンズオン

河合　洋介（国立国際医療研究センター）

2023年12月07日作成

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

[bcftools](https://github.com/samtools/bcftools)

[PSMC](https://github.com/lh3/psmc)

今回はスパコンにダウンロード済みのファイルを使います。

#### プログラムの確認
```
bcftools
psmc
```

ヘルプメッセージが表示されることを確認

### 使用するデータ
シミュレーションで作った*Fictus yaponesiae*のゲノムデータ(`FK01.bam`)を使います。

#### データの確認
* BAMファイル
```
samtools tview /home/bioarchaeology-pg/kawai/FK01.bam
```
表示がうまくいかない場合は
```
samtools view /home/bioarchaeology-pg/kawai/FK01.bam | less -S
```

* リファレンスファイル
```
less /home/bioarchaeology-pg/kawai/yaponesia_reference.fasta
```

### 解析ディレクトリの作成
```
# ホームディレクトに移動して作業ディレクトリを作成
cd
mkdir psmc
cd psmc
```

## 実行
### データの準備
* `bcftools`でBAMファイルからバリアント情報を抽出して、コンセンサス配列をfastq形式で抽出する。
```
bcftools mpileup -Ou -f /home/bioarchaeology-pg/kawai/yaponesia_reference.fasta /home/bioarchaeology-pg/kawai/FK01.bam | \
  bcftools call -c - | \
  vcfutils.pl vcf2fq -d 10 -D 100 | \
  gzip -c > FK01.fq.gz
```

>[!NOTE]
>コマンド入力でバックスラッシュ(`\`)に続けてエンターを入力すると改行される。

* 出力ファイルの確認
```
zless FK01.fq.gz
```
>[!NOTE]
>gzip圧縮されたデータは`zless`コマンドで中身を確認できる。


* PSMCに付属するツール(`fq2psmcfa`)で入力形式に変換する
```
fq2psmcfa FK01.fq.gz > FK01.psmcfa
```
* 出力ファイルの確認
```
less FK01.psmcfa
```
* PSMCの実行
```
psmc -t10 -p "10+5*3+4" -o FK01.psmc FK01.psmcfa
```
>[!NOTE]
>-t : 最大の世代数（2N<sub>0</sub>)
>

