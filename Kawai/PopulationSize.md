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
>-t : 最大の世代数（集団サイズ2N<sub>0</sub>が単位)
>
>-p : 有効集団サイズを推定する時間の分割パターンの指定。"10+5*3+4"は時間を10,3,3,3,3,3,4に分割して、それぞれの区画の有効手段サイズを推定する。
>
>-o : 出力ファイル名（.psmcをつけておく）

* 出力結果の確認
```
less -S FK01.psmc
```

* プロットの出力
```
psmc_plot.pl -u 2.0e-07 -g 1 -x 100 -L -p FK01 FK01.psmc
```
>[!NOTE]
>-u : 突然変異率
>
>-g : 世代時間（年）
>
>-x : 表示する最小の世代
>
>-p : PDFで出力（デフォルトはeps形式)

* 結果の確認

ターミナルで画像ファイルを閲覧することはできないので、ローカルの端末にPDFをダウンロードする。

ローカルのターミナルでscpコマンドで転送する方法
```
scp username-pg@gwa.ddbj.nig.ac.jp:/home/username-pg/psmc/FK01.pdf .
```
WinSCP(windows)やcyberduck(mac)で転送することも可能です。

## 実習
[ヒトゲノムの公共データ](https://sc.ddbj.nig.ac.jp/advanced_guides/advanced_guide_2023#ヒト全ゲノム解析の公共データの再解析データセット)を使ってPSMCで解析する。

遺伝研スパコンの中からは直接アクセス可能です。

`ls /usr/local/shared_data/public-human-genomes/GRCh38/`

* CRAMファイルからPSMCの入力ファイルを作成（時間がないのでハンズオンでは省略）
```
bcftools mpileup -Ou -f Homo_sapiens_assembly38.fasta \
  /usr/local/shared_data/public-human-genomes/GRCh38/1000Genomes/CRAM/NA12878/NA12878.cram | \
  bcftools call -c - | vcfutils.pl vcf2fq -d 10 -D 100 | gzip -c > CEU_NA12878.fq.gz
fq2psmcfa CEU_NA12878.fq.gz > CEU_NA12878.psmcfa
```
>[!NOTE]
>CRAMは圧縮率の高いマッピングファイルのフォーマットです。
>
>BAMの半分程度の大きさですが、開く時にリファレンス配列を参照する必要があります。

作成済みのpsmcfaファイルが`/home/bioarchaeology-pg/kawai/`に置いてあります。
```
CEU_NA12878.psmcfa
YRI_NA18496.psmcfa
JPT_NA18945.psmcfa
Karitiana_SS6004476.psmcfa
```



